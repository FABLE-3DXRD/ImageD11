#!/usr/bin/python
#!C:\Python24\python.exe
"""
Docstring
# By  Jaroslaw Butanowicz
# summer 2006
# modified by JPW from Activestate python recipies
#  - limited to localhost access for security
# Thanks also to AG..
"""

from SimpleXMLRPCServer import SimpleXMLRPCServer
import sys, logging, threading

class Server(SimpleXMLRPCServer):
    """
    (roughly) Recipes:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/165375
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/496700
    """
    accessList=('127.0.0.1')

    def __init__(self, *args, **kwds):
        SimpleXMLRPCServer.__init__(self, *args, **kwds)
        self.quit = False
        
    def server_bind(self):
        """
        Could be used to reuse same port (ouch)
        """
        # see unix $ man 2 setsockopt
        # Don't want to rebind I think - prefer to crash out if running
        # self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        SimpleXMLRPCServer.server_bind(self)

    def verify_request(self, request, client_address):
        """
        Not quite security...
        """
        if client_address[0] in self.accessList:
            return 1
        else:
            logging.debug("refused request from %s "% (
                           str(client_address)))
            return 0
        
#    def do_POST(self):
#        """
#        Removed - dont believe this is ever called by XMLRPC server
#        """
        
    def serve_forever(self):
        """ Mainloop of the server """
        logging.debug("Server serve_forever called")
        self.quit = False
        while not self.quit:
            self.handle_request()
            # try to help out the reads of stderr/stdout in the java gui
            sys.stdout.flush()
            sys.stderr.flush()


class XMLRPCThread(threading.Thread):
    """ Thread to run the server while we wait for control-c """
    def __init__(self, port):
        """ 
        Some things appear to be thread local in this namespace here
        This function will run the server in thread
        """            
        threading.Thread.__init__(self)
        self.port = port
        
    def run(self):
        """ The actual bit to spawn - server mainloop at end"""
        server = Server( ("localhost", self.port ) , logRequests = 0 )
        def shutdown():
            """ stop the server """
            logging.info( "ImageD11Server(): called shutdown()" ) 
            server.quit = True
        server.register_function(shutdown)
                
        def setLevel(level):
            """ Set the logging level of the root console logger"""
            logging.info( "ImageD11Server(): called setLevel "+str(level) ) 
            rt = logging.getLogger('')
            rt.setLevel(level)
            
        server.register_function(setLevel)
    
        from ImageD11 import guicommand
        guicommander = guicommand.guicommand()            
        server.register_instance(guicommander)
    
        #Go into the main listener loop
        logging.info("Python Server is listening on port " + str(self.port))
        server.serve_forever()
    
    

def RunXMLRPCServer(port):
    """Run the XMLRPCServer, but do so in a different thread.
    The main thread simply sleeps so that it can respond to Ctrl+Break
    http://mail.python.org/pipermail/python-list/2003-July/212751.html
    """
    
    # real function starts here
    import signal, time
    # the lambda binds func(port) to func() having no args
    # define a script/module level Handler which writes DEBUG messages or
    # higher to the sys.stdout
    th = XMLRPCThread(port)
    th.setDaemon(1)
    th.start()
    # Make ctrl+break raise a KeyboardInterrupt exception.
    if sys.platform == "win32":
        signal.signal(signal.SIGBREAK, signal.default_int_handler)
    else:
        logging.debug("signal.SIGBREAK not needed on this platform")
    try:
        while th.isAlive():
            time.sleep(1)
    except KeyboardInterrupt:
        logging.info( "Break pressed." )
        sys.stdout.flush()
        

if __name__=="__main__":

# drop everything on stdin/stdout
#    hdlr = logging.FileHandler('/tmp/xmlrpcserver.log')
#    formatter = logging.Formatter("%(asctime)s  %(levelname)s  %(message)s")
#    hdlr.setFormatter(formatter)
#    logger.addHandler(hdlr)


    use_this_port = 51234 # high randomish private number
    if len(sys.argv)>1:
        try:
            use_this_port = int(sys.argv[1])
        except (ValueError, IndexError) :
            print "usage: %s [portnumber]"% (sys.argv[0])
            
            
    # Set up the logging stuff
    console = logging.StreamHandler(sys.stdout)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-8s : %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    console.setLevel(logging.DEBUG)
    
    root = logging.getLogger('')
    root.addHandler(console)
    root.setLevel(logging.DEBUG) # should we process everything...?

    RunXMLRPCServer(use_this_port)



