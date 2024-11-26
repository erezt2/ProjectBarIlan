import scipy.io
import datetime, shutil
import secrets, itertools, random
import numpy as np
import time, os
from helpers import *
import chipwhisperer as cw


def posture(target, cmd, val, rd=None):
    # if simulate:
    #     return
    target.flush()
    target.simpleserial_write(cmd, val)
    if rd != None:
        return target.simpleserial_read('r', rd)

def wosture(target, cmd, val, rd):
    return get_bits(posture(target, cmd, get_bytes(val), rd))

def roll_seed():
    posture(target, 's', bytearray([random.randrange(256) for _ in range(4)]), 4)


def MEM(mem, **kwargs):
    return PUF(mem, None, roll_option="no-error", **kwargs)

def date_str():
    return datetime.datetime.now().strftime('%Y-%m-%d %H-%M')

def generate_file(env, g, *names):
    s = "\n"
    # g = globals()
    for nm in names:
        s += nm
        s += " = \n"
        s += str(g[nm])
        s += "\n"

    date_time = date_str()
    with open(f"{env.root}/{env.subfolder}/results/{env.name}.txt", 'w') as f:
        f.write(env.readme+"\n\n"+s)

# def traces(index):
#     trcs = []
#     for t in range(TRACES_NUM):
#         target.flush()
#         scope.arm()

#         msg = execute(index, t)
#         target.simpleserial_write('p', msg)
#         ret = scope.capture()
#         if ret:
#             print("ERROR")
#         traces = scope.get_last_trace()
#         fdb = target.simpleserial_read('r', 16)
#         ending(fdb, t)
#         trcs.append(traces)

#         # print(get_bytes(k)[0:8] == fdb[0:8]) 
#     return trcs

def write_append(file, data, append):
    if not append:
        scipy.io.savemat(file, {'traces': data})
    else:
        dat = np.concatenate((scipy.io.loadmat(file)['traces'], data), 0)
        scipy.io.savemat(file, {'traces': dat})

def traces(name, scope, target, env):
    trcs = []
    # for t in :
    #     for j, op in :
    a = range(env.TRACES_NUM)
    b = range(len(env.operations))
    l = list(itertools.product(a, b))
    random.shuffle(l)
    labels = []
    target.flush()
    for t, op in l:
        scope.adc.offset = 0
        msg = env.execute(env.index, env.operations[op], t)
        
        if env.flush:
            target.flush()
        traces = []
        for seg in range(env.segments):
            scope.arm()
            target.simpleserial_write('p', msg)
            ret = scope.capture()
            if ret:
                print(f"ERROR ON: {op = }, {t = }, {len(labels) = }")
                print(f"message: {msg}")
            traces = traces + list(scope.get_last_trace())
            if env.read > 0:
                fdb = target.simpleserial_read('r', env.read)
            else:
                fdb = None
            scope.adc.offset += scope.adc.samples
        
        throw = env.ending(fdb, env.operations[op], t)
        if not throw:
            trcs.append(traces)
            labels.append([op])
    file_n = f"{env.root}/{env.subfolder}/traces/{name}.mat"
    # temp_n = f"{env.root}/{env.subfolder}/traces/temporary.mat"
    if os.path.isfile(file_n):
        os.remove(file_n)
    scipy.io.savemat(file_n, {"traces": trcs, "labels": labels})
        # print(get_bytes(k)[0:8] == fdb[0:8]) 
    # os.rename(temp_n, file_n)
    print(name, end=" ")
    scope.adc.offset = 0
    return 0

def create_tasker(env, names, labels, out_name):
    folder = f"{env.root}/{env.subfolder}/traces"
    file_n = f"{folder}/{out_name}.mat"
    try:
        os.remove(file_n)
    except FileNotFoundError:
        pass
    scipy.io.savemat(file_n,  {"files": [f"{folder}/{nm}.mat" for nm in names], "labels":labels, "desc": env.desc, "title": env.title, "name": env.name, "graph": env.graph}) # the size is only correct if there is a constant amount of traces per set! change this line if this fact is changed

def split_file(env, names, labels, out_name):
    raise "Err"
    # remember to change env. desc + title + name before runs!
    # from the joint matrix of the files in names, bring out the labels in labels and merge into out_name
    folder = f"{env.root}/{env.subfolder}/traces"
    inv = {l: idx for idx, l in enumerate(labels)}
    n = len(labels)
    o_contents = [] # [] for i in range(n)
    o_labels = [] # [] for i in range(n)
    for name in names:
        dat = scipy.io.loadmat(f"{folder}/{name}.mat")
        lbl = dat["labels"]
        trc = dat["traces"]
        for i in range(len(lbl)):
            idx = inv.get(lbl[i][0])
            if idx is not None:
                o_contents.append(trc[i]) # [idx]
                o_labels.append([idx]) # [idx]
        del dat
        
    #for l in range(n):
    file_n = f"{folder}/{out_name}.mat" # operations[labels[l]]
    # temp_n = f"{folder}/temporary.mat" # operations[labels[l]]
    try:
        os.remove(file_n)
    except FileNotFoundError:
        pass
    scipy.io.savemat(file_n,  {"traces": o_contents, "labels": o_labels, "desc": env.desc, "title": env.title, "name": env.name, "graph": env.graph}) # [l]
    # os.rename(temp_n, file_n)

def change_graph(env, names, graph):
    folder = f"{env.root}/{env.subfolder}/traces"
    for name in names: # god FORBID there was some way to edit this file besides loading it and copying it
        file_n = f"{folder}/{name}.mat"
        dat = scipy.io.loadmat(file_n)
        os.remove(file_n)
        scipy.io.savemat(file_n,  {"traces": dat["traces"], "labels": dat["labels"], "desc": dat["desc"], "title": dat["title"], "name": dat["name"], "graph": graph})

def save_readme(env):
    folder = f"{env.root}/{env.subfolder}/results"
    file = f"{folder}/{env.name}.txt"
    with open(file, 'w') as f:
        f.write(env.readme)
    

def remove_files(env, names):
    folder = f"{env.root}/{env.subfolder}/traces/"
    for name in names:
        os.remove(f"{folder}/{name}.mat")
    


def load(scope, file, env):
    a = f"../hexes/{file}.hex"
    if env is not None:
        folder = f"{env.root}/{env.subfolder}/results"
        
        c = f"../hexes/{file}.c"
        p = f"./{env.jupyter}.ipynb"
        shutil.copy(c, folder)
        try:
            os.remove(f"{folder}/{env.name}.c")
        except FileNotFoundError:
            pass
        os.rename(f"{folder}/{file}.c", f"{folder}/{env.name}.c")
    
        shutil.copy(p, folder)
        try:
            os.remove(f"{folder}/{env.name}.ipynb")
        except FileNotFoundError:
            pass
        os.rename(f"{folder}/{env.jupyter}.ipynb", f"{folder}/{env.name}.ipynb")
    
    cw.program_target(scope, cw.programmers.STM32FProgrammer, a)
    
    


class ENV:
    def __init__(self):
        pass