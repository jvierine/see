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
                     "overlap_fraction":"2",                     
                     "ch":'["chc"]',
                     "offset":"1530",
                     "fscale":'"Hz"',
                     "fast":"false",
                     "n_avg":'"1"',
                     "debug":"false",
                     "debug_timing":"false",
                     "show_plot":"true",
                     "plot_carrier_pwr":"false",                     
                     "fmin":"-5e3",
                     "prefix":'"wb"',
                     "nsubsteps":"1",
                     "vmin":"-3",
                     "vmax":"60.0",
                     "n_cycles":"1",
                     "cycle_len":"1800.0",                     
                     "trim_end":"0",
                     "xc":"false",
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
        self.cycle_len=float(json.loads(c["config"]["cycle_len"]))        
        self.nsteps=int(json.loads(c["config"]["nsteps"]))
        self.nfft=int(json.loads(c["config"]["nfft"]))
        self.nsubsteps=int(json.loads(c["config"]["nsubsteps"]))        
        self.fmin=float(json.loads(c["config"]["fmin"]))
        self.fmax=float(json.loads(c["config"]["fmax"]))
        self.ch=json.loads(c["config"]["ch"])
        self.offset=int(json.loads(c["config"]["offset"]))
        self.trim_end=int(json.loads(c["config"]["trim_end"]))        
        self.f0=float(json.loads(c["config"]["f0"]))
        self.vmin=float(json.loads(c["config"]["vmin"]))
        self.vmax=float(json.loads(c["config"]["vmax"]))        
        self.fstep=float(json.loads(c["config"]["fstep"]))
        self.fscale=json.loads(c["config"]["fscale"])
        self.fast=bool(json.loads(c["config"]["fast"]))
        self.xc=bool(json.loads(c["config"]["xc"]))        
        self.debug=bool(json.loads(c["config"]["debug"]))
        self.plot_carrier_pwr=bool(json.loads(c["config"]["plot_carrier_pwr"]))        
        self.show_plot=bool(json.loads(c["config"]["show_plot"]))
        self.n_avg=int(json.loads(c["config"]["n_avg"]))
        self.overlap_fraction=float(json.loads(c["config"]["overlap_fraction"]))        
        self.n_cycles=int(json.loads(c["config"]["n_cycles"]))
        self.prefix=json.loads(c["config"]["prefix"])
        self.debug_timing=bool(json.loads(c["config"]["debug_timing"]))                

    def __str__(self):
        out="Configuration\n"
        for e in dir(self):
            if not callable(getattr(self,e)) and not e.startswith("__"):
                out+="%s = %s\n"%(e,getattr(self,e))
        return(out)
