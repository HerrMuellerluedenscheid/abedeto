#***************************************************************************** 
# client.py
#
# ArcLink client library
#
# (c) 2004 Andres Heinloo, GFZ Potsdam
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2, or (at your option) any later
# version. For more information, see http://www.gnu.org/
#*****************************************************************************

import bz2
import socket
from seiscomp import logs
from seiscomp.xmlparser import _MyContentHandler
from seiscomp.xmlparser import *

LINESIZE      = 1024
BLOCKSIZE     = 512

STATUS_UNSET  = 0
STATUS_PROC   = 1
STATUS_CANCEL = 2
STATUS_OK     = 3
STATUS_WARN   = 4
STATUS_ERROR  = 5
STATUS_RETRY  = 6
STATUS_DENIED = 7
STATUS_NODATA = 8

def arclink_status_string(s):
    if   s == STATUS_PROC:   return "PROCESSING"
    elif s == STATUS_CANCEL: return "CANCELLED"
    elif s == STATUS_OK:     return "OK"
    elif s == STATUS_WARN:   return "WARNING"
    elif s == STATUS_ERROR:  return "ERROR"
    elif s == STATUS_RETRY:  return "RETRY"
    elif s == STATUS_DENIED: return "DENIED"
    elif s == STATUS_NODATA: return "NODATA"
    elif s == STATUS_UNSET:  return "UNSET"

    return "UNKNOWN"

class _StatusAttr(XAttribute):
    def __init__(self, default = STATUS_UNSET):
        self._default = default

    def toxml(cls, val):
        return arclink_status_string(val)

    toxml = classmethod(toxml)

    def fromxml(self, val = None):
        if   val == "PROCESSING":  return STATUS_PROC
        elif val == "CANCELLED":   return STATUS_CANCEL
        elif val == "OK":          return STATUS_OK
        elif val == "WARNING":     return STATUS_WARN
        elif val == "ERROR":       return STATUS_ERROR
        elif val == "RETRY":       return STATUS_RETRY
        elif val == "DENIED":      return STATUS_DENIED
        elif val == "NODATA":      return STATUS_NODATA
        elif val == "UNSET":       return STATUS_UNSET
        elif val == None:          return STATUS_UNSET
        
        raise ValueError, "invalid status value: " + val

class _RequestLineElement(XElement):
    def __init__(self):
        XElement.__init__(self, { "content"     : StringAttr(),
                                  "status"      : _StatusAttr(),
                                  "size"        : IntAttr(),
                                  "message"     : StringAttr() })

class _RequestVolumeElement(XElement):
    def __init__(self):
        XElement.__init__(self, { "id"          : StringAttr(),
                                  "status"      : _StatusAttr(),
                                  "size"        : IntAttr(),
                                  "message"     : StringAttr(),
                                  "line"        : _RequestLineElement })

class _RequestElement(XElement):
    def __init__(self):
        XElement.__init__(self, { "id"          : StringAttr(),
                                  "user"        : StringAttr(),
                                  "institution" : StringAttr(),
                                  "label"       : StringAttr(),
                                  "type"        : StringAttr(),
                                  "args"        : StringAttr(),
                                  "size"        : IntAttr(),
                                  "message"     : StringAttr(),
                                  "ready"       : BoolAttr(),
                                  "error"       : BoolAttr(),
                                  "volume"      : _RequestVolumeElement })

class _ArclinkElement(XElement):
# FIXME: server does not add namespace yet
#   xmlns = "http://geofon.gfz-potsdam.de/ns/arclink/0.1"
    def __init__(self):
        XElement.__init__(self, { "request"     : _RequestElement })

class ArclinkError(Exception):
    pass

class ArclinkCommandNotAccepted(ArclinkError):
    pass

class ArclinkAuthFailed(ArclinkError):
    def __init__(self, dcname):
        self.dcname = dcname
        ArclinkError.__init__(self, "authentication failed")
    
class ArclinkXMLError(ArclinkError):
    pass

class ArclinkTimeout(ArclinkError):
    pass

class Arclink(object):
    def __init__(self):
        self.__sock = None
        self.__fd = None
        self.__errstate = False
    
    def send_command(self, s, check_ok=True):
        self.__fd.write(s + "\r\n")
        self.__fd.flush()
        if not check_ok:
            return

        r = self.__fd.readline(LINESIZE).rstrip()
        if r == "OK":
            return
        elif r == "ERROR":
            raise ArclinkCommandNotAccepted, "command not accepted"
        else:
            raise ArclinkError, "unexpected response: " + r

    def get_errmsg(self):
        self.send_command("SHOWERR", False)
        return self.__fd.readline(LINESIZE).rstrip()

    def open_connection(self, host, port, user, passwd = None, inst = None,
        timeout = None, user_ip = None):
        self.__sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.__sock.settimeout(timeout)
        
        self.__sock.connect((host, port))
        self.__fd = self.__sock.makefile('rb+')
        try:
            self.send_command("HELLO", False)
            r = self.__fd.readline(LINESIZE).rstrip()
            if r == "ERROR":
                raise ArclinkCommandNotAccepted, "command not accepted"

            self.software = r
            self.organization = self.__fd.readline(LINESIZE).rstrip()
        
            if user == None:
                raise ArclinkAuthFailed, self.organization
        
            try:
                if(passwd == None):
                    self.send_command("USER " + user)
                else:
                    self.send_command("USER " + user + " " + passwd)
            except ArclinkCommandNotAccepted:
                raise ArclinkAuthFailed, self.organization

            if inst != None:
                self.send_command("INSTITUTION " + inst)

            if user_ip != None:
                self.send_command("USER_IP " + user_ip)

        except:
            self.__fd.close()
            self.__sock.close()
            raise

    def close_connection(self):
        self.__fd.close()
        self.__sock.close()
    
    def get_status(self, req_id = "ALL"):
        self.send_command("STATUS " + req_id, False)
        p = xml.sax.make_parser()
        element_stack = [XElement({ "arclink" : _ArclinkElement })]
        p.setFeature(xml.sax.handler.feature_namespaces, True)
        p.setContentHandler(_MyContentHandler(p, element_stack, True))
 
        try:
            line = self.__fd.readline(LINESIZE).rstrip()
            while line != "END":
                if not line:
                    raise ArclinkError, "unexpected end of XML data"
                elif line == "ERROR":
                    raise ArclinkError, self.get_errmsg()

                p.feed(line)
                line = self.__fd.readline(LINESIZE).rstrip()

            p.close()

        except xml.sax.SAXException, e:
            raise ArclinkXMLError, str(e)
        
        if(len(element_stack[0].arclink) != 1):
            raise ArclinkXMLError, "bad XML root element"
    
        return element_stack[0].arclink[0]

    def submit(self, reqf):
        line = reqf.readline(LINESIZE).strip()
        while not line:
            line = reqf.readline(LINESIZE).strip()
        
        if line.split()[0].upper() == "LABEL":
            self.send_command(line)
            line = reqf.readline(LINESIZE).strip()
        
        if line.split()[0].upper() != "REQUEST":
            raise ArclinkError, "incorrectly formulated request"
            
        self.send_command(line)
        
        line = reqf.readline(LINESIZE).strip()
        while line and line.split()[0].upper() != "END":
            self.send_command(line, False)
            line = reqf.readline(LINESIZE).strip()

        if len(line.split()) != 1:
            raise ArclinkError, "incorrectly formulated request"
        
        line = reqf.readline(LINESIZE)
        while line:
            if line.strip():
                raise ArclinkError, "incorrectly formulated request"

            line = reqf.readline(LINESIZE)

        self.send_command("END", False)
        r = self.__fd.readline(LINESIZE).rstrip()
        if r == "ERROR":
            raise ArclinkError, self.get_errmsg()

        return r

    def __init_download(self, req_id, vol_id, pos, timeout):
        if vol_id == None:
            req_vol = req_id
        else:
            req_vol = req_id + "." + vol_id
        
        if pos == None:
            pos_ext = ""
        else:
            pos_ext = " " + pos
        
        saved_tmo = self.__sock.gettimeout()
        try:
            if timeout == 0:
                cmdstr = "DOWNLOAD "
            else:
                cmdstr = "BDOWNLOAD "
                self.__sock.settimeout(timeout)
            
            self.send_command(cmdstr + req_vol + pos_ext, False)
            
            try:
                r = self.__fd.readline(LINESIZE).rstrip()
                if r == "ERROR":
                    raise ArclinkError, self.get_errmsg()
            
            except socket.timeout:
                if timeout == 0:
                    raise
                
                raise ArclinkTimeout, "timeout"

        finally:
            self.__sock.settimeout(saved_tmo)
            
        try:
            size = int(r)
        except ValueError:
            raise ArclinkError, "unexpected response: " + r

        return size
    
    def download_data(self, outfd, req_id, vol_id = None, pos = None,
        timeout = 0, decomp = None):
        
        bytes_read = 0
        size = self.__init_download(req_id, vol_id, pos, timeout)
        self.__errstate = False

        try:
            if decomp is None or decomp == "none":
                while bytes_read < size:
                    buf = self.__fd.read(min(BLOCKSIZE, size - bytes_read))
                    outfd.write(buf)
                    bytes_read += len(buf)
            elif decomp == "bzip2":
                z = bz2.BZ2Decompressor()
        
                while bytes_read < size:
                    zbuf = self.__fd.read(min(BLOCKSIZE, size - bytes_read))
                    buf = z.decompress(zbuf)
                    outfd.write(buf)
                    bytes_read += len(zbuf)
            else:
                raise ArclinkError, "unsupported compression"

            r = self.__fd.readline(LINESIZE).rstrip()
            if r != "END":
                raise ArclinkError, "END not found"

        finally:
            if bytes_read and bytes_read != size:
                self.__errstate = True

    def download_xml(self, db, req_id, vol_id = None, pos = None,
        timeout = 0, decomp = None):

        p = db.make_parser()

        bytes_read = 0
        size = self.__init_download(req_id, vol_id, pos, timeout)
        self.__errstate = False

        try:
            try:
                if decomp is None or decomp == "none":
                    while bytes_read < size:
                        buf = self.__fd.read(min(BLOCKSIZE, size - bytes_read))
                        p.feed(buf)
                        bytes_read += len(buf)

                elif decomp == "bzip2":
                    z = bz2.BZ2Decompressor()
        
                    while bytes_read < size:
                        zbuf = self.__fd.read(min(BLOCKSIZE, size - bytes_read))
                        buf = z.decompress(zbuf)
                        p.feed(buf)
                        bytes_read += len(zbuf)
                else:
                    raise ArclinkError, "unsupported compression"

                p.close()

            except xml.sax.SAXException, e:  # FIXME: not applicable to ElementTree
                raise ArclinkXMLError, str(e)
            
            r = self.__fd.readline(LINESIZE).rstrip()
            if r != "END":
                raise ArclinkError, "END not found"

        finally:
            if bytes_read and bytes_read != size:
                self.__errstate = True

    def errstate(self):
        return self.__errstate
    
    def purge(self, req_id):
        try:
            self.send_command("PURGE " + req_id)
        except ArclinkCommandNotAccepted:
            raise ArclinkError, self.get_errmsg()

