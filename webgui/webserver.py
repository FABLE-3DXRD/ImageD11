#Copyright Jon Berg , turtlemeat.com

import string,cgi,time
from os import curdir, sep
from BaseHTTPServer import BaseHTTPRequestHandler, HTTPServer
#import pri



class worker(object):

    def __init__(self):
        self.stuff = {"user":"Jon"}

    def makemenu(self):
        return """<html>
<body>
 Hello Matey
 <p>
 <form action="setStuff" method="get">
   User Name: <input name="user" type="text" />
   Currently %s
   <input type="submit" />
 </form>
 <form action="run" method="get">
   OK, do it, do it now
   <input type="submit" />
 </form>
</body>
</html>
    """%(self.stuff["user"])

    def set(self, k, v):
        self.stuff[k] = v



myworker = worker()


class MyHandler(BaseHTTPRequestHandler):

    def do_GET(self):
        print "Got a GET"
        if "setStuff" in self.path:
            items = [ x.split("=") for x in self.path.split("?")[1:] ]
            print items
            for k,v in items:
                myworker.set(k,v)
        try:
            self.send_response(200)
            self.send_header('Content-type',	'text/html')
            self.end_headers()
            if self.path.find("favicon.ico")>-1:
                self.wfile.write( open("favicon.ico").read() )
            else:
                self.wfile.write(myworker.makemenu())

            print "hi", self.path
            return
                
        except IOError:
            self.send_error(404,'File Not Found: %s' % self.path)
     

    def do_POST(self):
        global rootnode
        print "Got a POST"
        try:
            ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))
            if ctype == 'multipart/form-data':
                query=cgi.parse_multipart(self.rfile, pdict)
            self.send_response(301)
            
            self.end_headers()
            upfilecontent = query.get('upfile')
#            print "filecontent", upfilecontent[0]
            self.wfile.write("<HTML>POST OK.<BR><BR>");
            self.wfile.write(upfilecontent[0]);
            
        except :
            pass

def main():
    try:
        server = HTTPServer(('', 80), MyHandler)
        print dir(server)
        print server.server_port
        print 'started httpserver...'
        server.serve_forever()
    except KeyboardInterrupt:
        print '^C received, shutting down server'
        server.socket.close()

if __name__ == '__main__':
    main()

