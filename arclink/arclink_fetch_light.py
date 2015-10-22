#*****************************************************************************
# arclink_fetch.py
#
# ArcLink command-line client with routing support
#
# (c) 2009 Andres Heinloo, GFZ Potsdam
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2, or (at your option) any later
# version. For more information, see http://www.gnu.org/
#*****************************************************************************

import os
import sys
import datetime
import shutil
from optparse import OptionParser
from tempfile import TemporaryFile
from breqfast import BreqParser
from seiscomp import logs
from seiscomp.arclink.manager import *
from seiscomp.db import DBError
from seiscomp.db.generic.inventory import Inventory
from seiscomp.mseedlite import Input as MSeedInput, MSeedError
from seiscomp.fseed import SEEDVolume, SEEDError, _WaveformData

VERSION = "1.0 (2011.136)"

ORGANIZATION = "WebDC"
LABEL = "WebDC SEED Volume"

verbosity = 1

class SeedOutput(object):
    def __init__(self, fd, inv, resp_dict):
        self.__fd = fd
        self.__inv = inv
        self.__resp_dict = resp_dict
        self.__mseed_fd = TemporaryFile()

    def write(self, data):
        self.__mseed_fd.write(data)

    def close(self):
        try:
            try:
                seed_volume = SEEDVolume(self.__inv, ORGANIZATION, LABEL,
                    self.__resp_dict)

                self.__mseed_fd.seek(0)
                for rec in MSeedInput(self.__mseed_fd):
                    seed_volume.add_data(rec)

                seed_volume.output(self.__fd)

            except (MSeedError, SEEDError, DBError), e:
                logs.error("error creating SEED volume: " + str(e))

        finally:
            self.__mseed_fd.close()
            self.__fd.close()

class MSeed4KOutput(object):
    def __init__(self, fd):
        self.__fd = fd
        self.__mseed_fd = TemporaryFile()

    def write(self, data):
        self.__mseed_fd.write(data)

    def close(self):
        try:
            try:
                wfd = _WaveformData()

                self.__mseed_fd.seek(0)
                for rec in MSeedInput(self.__mseed_fd):
                    wfd.add_data(rec)

                wfd.output_data(self.__fd, 0)

            except (MSeedError, SEEDError, DBError), e:
                logs.error("error reblocking Mini-SEED data: " + str(e))

        finally:
            self.__mseed_fd.close()
            self.__fd.close()

def show_status(rqstat):
    if rqstat.error:
        req_status = "ERROR"
    elif rqstat.ready:
        req_status = "READY"
    else:
        req_status = "PROCESSING"

    logs.info("Request ID: %s, Label: %s, Type: %s, Args: %s" % \
        (rqstat.id, rqstat.label, rqstat.type, rqstat.args))
    logs.info("Status: %s, Size: %d, Info: %s" % \
        (req_status, rqstat.size, rqstat.message))

    for vol in rqstat.volume:
        logs.info("    Volume ID: %s, Status: %s, Size: %d, Info: %s" % \
            (vol.id, arclink_status_string(vol.status), vol.size, vol.message))

        for rqln in vol.line:
            logs.info("        Request: %s" % (rqln.content,))
            logs.info("        Status: %s, Size: %d, Info: %s" % \
              (arclink_status_string(rqln.status), rqln.size, rqln.message))

    logs.info("")

def parse_native(req, input_file):
    fd = open(input_file)
    try:
        rqline = fd.readline()
        while rqline:
            rqsplit = rqline.split()
            if len(rqsplit) < 3:
                logs.error("invalid request line: %s" % (rqline,))
                rqline = fd.readline()
                continue

            try:
                start_time = datetime.datetime(*map(int, rqsplit[0].split(",")))
                end_time = datetime.datetime(*map(int, rqsplit[1].split(",")))
            except ValueError, e:
                logs.error("syntax error (%s): %s" % (str(e), rqline))
                rqline = fd.readline()
                continue

            network = rqsplit[2]
            station = "."
            channel = "."
            location = "."

            i = 3
            if len(rqsplit) > 3 and rqsplit[3] != ".":
                station = rqsplit[3]
                i += 1
                if len(rqsplit) > 4 and rqsplit[4] != ".":
                    channel = rqsplit[4]
                    i += 1
                    if len(rqsplit) > 5 and rqsplit[5] != ".":
                        location = rqsplit[5]
                        i += 1
                        
            while len(rqsplit) > i and rqsplit[i] == ".":
                i += 1
            
            constraints = {}
            for arg in rqsplit[i:]:
                pv = arg.split('=', 1)
                if len(pv) != 2:
                    raise ArclinkHandlerError, "invalid request syntax"
                
                constraints[pv[0]] = pv[1]

            req.add(network, station, channel, location, start_time, end_time,
                constraints)
    
            rqline = fd.readline()

    finally:
        fd.close()

def parse_breqfast(req, input_file):
    parser = BreqParser()
    parser.parse_email(input_file)
    req.content = parser.reqlist
    if parser.failstr:
        logs.error(parser.failstr)

def add_verbosity(option, opt_str, value, parser):
    global verbosity
    verbosity += 1

def add_quietness(option, opt_str, value, parser):
    global verbosity
    verbosity -= 1

def process_options():
    parser = OptionParser(usage="usage: %prog [-a host:port] [-f format] [-k kind] [-n] [-g] [-t timeout] [-x retries] [-v] [-q] -u user -o file request",
      version="%prog v" + VERSION)

    parser.set_defaults(address = "webdc.eu:18001",
                        request_format = "native",
                        data_format = "mseed",
                        no_resp_dict = False,
                        rebuild_volume = False,
                        timeout = 300,
                        retries = 5)

    parser.add_option("-a", "--address", type="string", dest="address",
      help="address of primary ArcLink node (default %default)")

    parser.add_option("-f", "--request-format", type="string", dest="request_format",
      help="request format: breqfast, native (default %default)")

    parser.add_option("-k", "--data-format", type="string", dest="data_format",
      help="data format: mseed, mseed4k, fseed, dseed, inv[entory] (default %default)")

    parser.add_option("-n", "--no-resp-dict", action="store_true", dest="no_resp_dict",
      help="avoid using response dictionary (default %default)")

    parser.add_option("-g", "--rebuild-volume", action="store_true", dest="rebuild_volume",
      help="rebuild SEED volume (default %default)")

    parser.add_option("-t", "--timeout", type="int", dest="timeout",
      help="timeout in seconds (default %default)")

    parser.add_option("-x", "--retries", type="int", dest="retries",
      help="download retries (default %default)")

    parser.add_option("-v", action="callback", callback=add_verbosity,
      help="increase verbosity level")
    
    parser.add_option("-q", action="callback", callback=add_quietness,
      help="decrease verbosity level")
    
    parser.add_option("-u", "--user", type="string", dest="user",
      help="user's e-mail address")

    parser.add_option("-o", "--output-file", type="string", dest="output_file",
      help="file where downloaded data is written")

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error("incorrect number of arguments")

    if options.user == None:
        parser.error("username required")
    
    if options.output_file == None:
        parser.error("output file required")
    
    if options.data_format.upper() != "FSEED" and options.rebuild_volume:
        parser.error("-g is only applicable to FSEED format")
    
    return (options.address, options.request_format, options.data_format,
      not options.no_resp_dict, options.rebuild_volume, options.user,
      options.timeout, options.retries, options.output_file, args[0])

    parser.set_defaults(address = "webdc.eu:18001",
                        request_format = "native",
                        data_format = "mseed",
                        no_resp_dict = False,
                        rebuild_volume = False,
                        timeout = 300,
                        retries = 5)
def main(user, output_file, input_file, addr='webdc.eu:18001', request_format='native', data_format='mseed',
         resp_dict=False, rebuild_volume=False, timeout=300, retries=5):
    reblock_mseed = False
    use_inventory = False
    use_routing = True

    req_args = {"compression": "bzip2"}
    if data_format.upper() == "MSEED":
        req_type = "WAVEFORM"
        req_args["format"] = "MSEED"

    elif data_format.upper() == "MSEED4K":
        req_type = "WAVEFORM"
        req_args["format"] = "MSEED"
        reblock_mseed = True

    elif data_format.upper() == "FSEED":
        req_type = "WAVEFORM"
        if rebuild_volume:
            req_args["format"] = "MSEED"
        
        else:
            req_args["format"] = "FSEED"
    
    elif data_format.upper() == "DSEED":
        req_type = "RESPONSE"
        use_routing = False
    
    elif len(data_format) >= 3 and data_format.upper() == "INVENTORY"[:len(data_format)]:
        req_type = "INVENTORY"
        req_args["instruments"] = "true"
        use_routing = False

    else:
        logs.error("unsupported data format: %s" % (data_format,))
        return 1
    
    if resp_dict:
        req_args["resp_dict"] = "true"
    else:
        req_args["resp_dict"] = "false"
    
    fd_out = open(output_file, "wb")
    mgr = ArclinkManager(addr, user, socket_timeout=timeout, download_retry=retries)
    req = mgr.new_request(req_type, req_args)

    if request_format == "native":
        parse_native(req, input_file)

    elif request_format == "breqfast":
        parse_breqfast(req, input_file)

    else:
        logs.error("unsupported request format: %s" % (request_format,))
        return 1

    if not req.content:
        logs.error("empty request")
        return 1
    
    wildcards = False
    for i in req.content:
        for j in i[:4]:
            if j.find("*") >= 0 or j.find("?") >= 0:
                wildcards = True
                break

    if rebuild_volume or wildcards:
        use_inventory = True

    try:
        (files, inv, req_ok, req_noroute, req_nodata) = mgr.get_data(req, use_inventory, use_routing)

    except ArclinkError, e:
        logs.error(str(e))
        return

    if verbosity > 1:
        logs.info("the following data requests were sent:")
        for req in req_ok:
            logs.info(req.dcname)
            
            try:
                show_status(req.status())

            except ArclinkError, e:
                logs.error(str(e))

    if verbosity > 0:
        if req_nodata:
            logs.info("the following entries returned no data:")
            req_nodata.dump(sys.stdout)

        if req_noroute:
            logs.info("the following entries could not be routed:")
            req_noroute.dump(sys.stdout)

    if rebuild_volume:
        logs.info("rebuilding SEED volume")
        fd_out = SeedOutput(fd_out, inv, resp_dict)
    
    elif reblock_mseed:
        logs.info("reblocking Mini-SEED data")
        fd_out = MSeed4KOutput(fd_out)

    for fd in files:
        fd.seek(0)
        shutil.copyfileobj(fd, fd_out)
        fd.close()

    fd_out.close()

def _debug(s):
    if verbosity > 1:
        print s
        sys.stdout.flush()

def _info(s):
    if verbosity > 0:
        print s
        sys.stdout.flush()

if __name__ == "__main__":
    logs.info = _info
    logs.debug = _debug
    sys.exit(main())

