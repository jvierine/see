try:
    import configparser
except ImportError as e:
    import configparser2 as configparser
    
import json
import os 
class see_config:
    def __init__(self,fname=None):
        c=configparser.ConfigParser()
        # initialize with default values
        c["config"]={"t0":'"1602761400.0"',
                     "sample_rate":"1000000.0",
                     "center_freq":"4.04e6",
                     "data_dirs":'["/data0/2020.10.15/test1e6_4.04e6","/data1/2020.10.15/test1e6_4.04e6"]',
                     "step_len":"10",
                     "nsteps":"133", 
                     "fstep":"3.125e3",
                     "f0":"3.9e6",
                     "nfft":"262144",
                     "ch":'["chc","chd"]',
                     "offset":"1530",
                     "fmin":"-5e3",
                     "fmax":"5e3"}

        if fname != None:
            if os.path.exists(fname):
                print("reading %s"%(fname))
                c.read(fname)
            else:
                print("configuration file %s doesn't exist. using default values"%(fname))
        self.fname=fname
        self.t0=float(json.loads(c["config"]["t0"]))
        self.sample_rate=float(json.loads(c["config"]["sample_rate"]))
        self.center_freq=float(json.loads(c["config"]["center_freq"]))
        self.data_dirs=json.loads(c["config"]["data_dirs"])
        self.step_len=float(json.loads(c["config"]["step_len"]))
        self.nsteps=int(json.loads(c["config"]["nsteps"]))
        self.nfft=int(json.loads(c["config"]["nfft"]))
        self.fmin=float(json.loads(c["config"]["fmin"]))
        self.fmax=float(json.loads(c["config"]["fmax"]))
        self.ch=json.loads(c["config"]["ch"])
        self.offset=int(json.loads(c["config"]["offset"]))
        self.f0=float(json.loads(c["config"]["f0"]))
        self.fstep=float(json.loads(c["config"]["fstep"]))

    def __str__(self):
        out="Configuration\n"
        for e in dir(self):
            if not callable(getattr(self,e)) and not e.startswith("__"):
                out+="%s = %s\n"%(e,getattr(self,e))
        return(out)