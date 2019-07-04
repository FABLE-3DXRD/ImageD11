
import sys, os


PYT = sys.executable
FILT= os.path.join("..","..","scripts","filtergrain.py")
FIT = os.path.join("..","..","scripts","fitgrain.py")



os.system("%s %s %s"%(PYT, FILT, 
        "-f peaks.out_merge_t100 -u g3.ubi -F new.flt -p g3.pars -g 0 -t 0.15"))

os.system("%s %s %s"%(PYT, FIT, 
        "-u g3.ubi -f new.flt -p g3.pars -P new.pars -U new.ubi -S 300"))

os.system("%s %s %s"%(PYT, FILT, 
        "-f peaks.out_merge_t100 -u new.ubi -F new.flt -p new.pars -g 0 -t 0.1"))

os.system("%s %s %s"%(PYT, FIT, 
        "-u g3.ubi -f new.flt -p new.pars -P new.pars -U new.ubi -S 300"))

os.system("%s %s %s"%(PYT, FILT, 
        "-f peaks.out_merge_t100 -u new.ubi -F new.flt -p new.pars -g 0 -t 0.05"))

os.system("%s %s %s"%(PYT, FIT, 
        "-u g3.ubi -f new.flt -p new.pars -P new.pars -U new.ubi -S 300"))

os.system("%s %s %s"%(PYT, FILT, 
        "-f peaks.out_merge_t100 -u new.ubi -F new.flt -p new.pars -g 0 -t 0.025"))

os.system("%s %s %s"%(PYT, FIT, 
        "-u new.ubi -f new.flt -p new.pars -P new.pars -U new.ubi -S 300"))

os.system("%s %s %s"%(PYT, FILT, 
        "-f peaks.out_merge_t100 -u new.ubi -F new.flt -p new.pars -g 0 -t 0.025"))

os.system("%s %s %s"%(PYT, FIT, 
        "-u new.ubi -f new.flt -p new.pars -P new.pars -U new.ubi -S 300"))

