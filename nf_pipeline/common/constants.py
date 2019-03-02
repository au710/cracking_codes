import yaml
import os

def read_config(fname, fallback=None):

    if os.path.exists(fname):
        source = fname
    else:
        source = fallback

    with open(source, 'r') as stream:
        content = yaml.load(stream)
        return content


config = read_config("constants.yml")

Kap_raw = config["Kap_RAW"]
Ga_raw  = config["Ga_RAW"]
T_raw   = config["T_RAW"]
K_raw   = config["K_RAW"]
G_raw   = config["G_RAW"]
alp_raw = config["alp_RAW"]
phi_raw = config["phi_RAW"]
TIM     = config["TIM"]
N_t     = config["N_t"]
sig_raw = config["sig_RAW"]
