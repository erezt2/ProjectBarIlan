import random, math
from functools import reduce

def app(b):
    return None

default_flip = False

#def multi_op(func):
#    def wrapper(*x):
#        z = x[0]
#        for i in range(1,len(x)):
#            z = func(z, x[i])
#        return z
#    return wrapper
#
#@multi_op
def xor_aux(x,y):
    z = []
    if len(x) != len(y):
        print(len(x), len(y))
        raise Exception("mimatching lengths")
    for i in range(len(x)):
        z.append(x[i]^y[i])
    return z

def rand_weight(size, weight):
    pz = [0]*(size-weight) + [1]*weight
    return random.sample(pz, size)

def error_on(lst, pos, set=None):
    lst = list(lst)
    if pos<0:
        return lst
    if set is None:
        lst[pos] ^= 1
        return lst
    lst[pos] = set
    return lst
        

def bitwise_not(x):
    return [1-b for b in x]

def xor(*args):
    return reduce(xor_aux, args)


def lshift(lst, sft):
    return lst[sft:] + lst[0:sft]

def rshift(lst, sft):
    return lst[-sft:] + lst[0:-sft]
    
def get_bits_aux(x, flip=default_flip):
    x = int(x)
    l = []
    for i in range(8):
        l.append(x & 1)
        x >>= 1
    if flip:
        l = l[::-1]
    return l
    
def get_bits(z, flip=default_flip):
    return reduce(lambda x,y: x+y,list(map(lambda x_: get_bits_aux(x_, flip=flip), z)),[])

def get_byte_aux(x, flip=default_flip):

    # if len(x)!=8: not needed, functions the same
    #    x = x + (8-len(x))*[0]
    #    print("hehehe")
    if not flip:
        x=x[::-1]
    res = 0
    for y in x:
        res <<= 1
        res += y
    return res

def get_bytes(z, flip=default_flip):
    n = math.ceil((len(z)/8))
    p = [z[8*i:8*i+8] for i in range(n)]
    return bytearray(map(lambda x: get_byte_aux(x, flip=flip),p))

class PUF: #length-127
    stable_err = [10**-6, 10**-5]
    unstable_err = [10**-3, 3*10**-2]

    def generate_stable_mask(self, per_unstable):
        res = []
        for _ in range(self.len): # 1-stable, 0-unstable
            if per_unstable is None:
                t = 0
            elif random.random() < per_unstable:
                t = random.uniform(*self.unstable_err)
                self.unstable += 1
            else:
                t = random.uniform(*self.stable_err)
            res.append(t)
        return res

    def regenerate_mask(self, per_unstable):
        self.unstable = 0
        self.fault = self.generate_stable_mask(per_unstable)

    def __init__(self, puf: list, per_unstable=None, /,roll_option="majority7", uni=False):
        # puf definition
        self.puf = puf
        self.len = len(puf)

        # puf sample options
        self.unstable = 0
        self.roll_option = roll_option
        self.fault = self.generate_stable_mask(per_unstable)
        
        # puf fault options
        self.error = 0
        self.unidirectional = uni
        self.puf_raise = False
        
        # puf meta-data
        self.error_flag = False
        

    def roll_bit(self, index):
        # if self.fault[i] == 0: # this line doesn't change anything, just for semantics
        #     return self.puf[i]
        
        i = index
        if self.roll_option == "majority7":
            self.roll_option = "once"
            res = sum([self.roll_bit(i) for _ in range(7)], 0)
            self.roll_option = "majority7"
            return int(res > 3)
        elif self.roll_option == "once":
            return self.puf[i] ^ (random.random() < self.fault[i])
        elif self.roll_option == "no-error":
            return self.puf[i]
        
        raise Exception("invalid roll option")


    def roll(self):
        msg = []
        self.error_flag = False
        for i in range(self.len):
            b = self.roll_bit(i)
            if b != self.puf[i]:
                self.error_flag = True
            msg.append(b)
        return msg
    

    ### could change according to fault wanted ###

    def fault_insertion(self, msg, data): 
        if data == None:
            return msg
        i = data["index"]
        if i <= -1:
            return msg
        if self.error > random.random():
            if self.puf_raise:
                msg[i] = 1
        else:
            if self.unidirectional:
                msg[i] ^= 1
            else:
                msg[i] = 0
        return msg

    def sample(self, data=None, byte=True):
        msg = self.roll()
        msg = self.fault_insertion(msg, data)
        if not byte:
            return list(msg)
        return get_bytes(msg)
    
    #############################################