#!/usr/bin/python2.4

import SimpleXMLRPCServer

import guicommand

if __name__=="__main__":
    #server object
    guicommander = guicommand.guicommand()

    #merger_object.readpeaks(testfile)
    #    object.harvestpeaks()
    #merger_object.mergepeaks()
    #merger_object.filter()
    #merger_object.savepeaks(fltfile)
    #    sys.exit(0)
    server = SimpleXMLRPCServer.SimpleXMLRPCServer(("localhost", 8888))
    #server.register_instance(sum_object)

    server.register_instance(guicommander)

    #Go into the main listener loop
    print "Listening on port 8888"
    server.serve_forever()
