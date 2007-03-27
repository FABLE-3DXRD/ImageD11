#!/usr/bin/python2.4

# By  Jaroslaw Butanowicz
# summer 2006
# modified by JPW from Activestate python recipies
#  - limited to localhost access for security
#  


from SimpleXMLRPCServer import SimpleXMLRPCServer
import logging 

class Server(SimpleXMLRPCServer):
    """
    (roughly) Recipes:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/165375
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/496700
    """
    accessList=('127.0.0.1')

    def __init__(self,*args,**kwds):
        SimpleXMLRPCServer.__init__(self,*args,**kwds)
        
    def server_bind(self):
        # see unix $ man 2 setsockopt
        # Don't want to rebind I think - prefer to crash out if running
        # self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        SimpleXMLRPCServer.server_bind(self)

    def verify_request(self,request, client_address):
        """
        Not quite security...
        """
        if client_address[0] in self.accessList:
            return 1
        else:
            return 0
        
    def do_POST(self):
        clientIP, port = self.client_address
	    # Log client IP and Port
        logger.info('Client IP: %s - Port: %s' % (clientIP, port))
        try:
            # get arguments
            data = self.rfile.read(int(self.headers["content-length"]))
            # Log client request
            logger.info('Client request: \n%s\n' % data)
        
            response = self.server._marshaled_dispatch(
                    data, getattr(self, '_dispatch', None)
                )
	        # Log server response
            logger.info('Server response: \n%s\n' % response)
        
        except: # This should only happen if the module is buggy
            # internal error, report as HTTP server error
            self.send_response(500)
            self.end_headers()
        else:
            # got a valid XML RPC response
            self.send_response(200)
            self.send_header("Content-type", "text/xml")
            self.send_header("Content-length", str(len(response)))
            self.end_headers()
            self.wfile.write(response)

            # shut down the connection
            self.wfile.flush()
            self.connection.shutdown(1)



if __name__=="__main__":

    logger = logging.getLogger('xmlrpcserver')
    hdlr = logging.FileHandler('xmlrpcserver.log')
    formatter = logging.Formatter("%(asctime)s  %(levelname)s  %(message)s")
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)

    from ImageD11 import guicommand
    guicommander = guicommand.guicommand()

    import sys
    port = 51234
    if len(sys.argv)>1:
        try:
            port = int(sys.argv[1])
        except:
            print "usage: %s [portnumber]"%(sys.argv[0])

    server = Server( ("localhost", port ), logRequests=1 )
    server.register_instance(guicommander)

    #Go into the main listener loop
    print "Listening on port ",port
    server.serve_forever()
