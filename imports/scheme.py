import bchlib, random # this module gives up when decoding an error too big with the syndromes.
from helpers import *

swap = True # should be true
class Scheme:
    def __init__(self, scheme=None):
        self.b = bchlib.BCH(10, m=7, swap_bits=swap)
        self.helper = {}
        self.meta = {}
        self.scheme = scheme
        self.info = []


    def scheme_1(self, r_):  # r_ = r + e
        m_ = [random.randrange(2) for _ in range(64)]
        c_ = self.encode(m_)
        y = xor(xor(self.helper["enc"], c_), r_) # helper = w
        y__ = self.decode(y)
        k_ = xor(y__,m_)
        return k_

    def scheme_2(self, r_): # r_ = r + e
        m = [random.randrange(2) for _ in range(64)]
        c = self.encode(m)
        y = xor(c, r_)
        s_ = self.syndrome(y)
        s = xor(s_, self.helper["syn"]) # helper = w
        e = self.syn2error(s)
        r = xor(e, r_)
        return r[63:]

    def transpose(self, arr2d):
        return [[arr2d[j][i] for j in range(len(arr2d))] for i in range(len(arr2d[0]))]

    def xor2d(self, x, y):
        if len(x) != len(y):
            raise Exception("mimatching lengths")
        return [xor(x[i],y[i]) for i in range(len(x))]
    
    def repetition7_ste(self, syn):
        if len(syn) != 6:
            raise Exception("bad syndrome")
        if sum(syn) > 3:
            return [1] + [x^1 for x in syn]
        else:
            return [0] + [x^0 for x in syn]

    def stacked_r7ste(self, syn7):
        xx = self.transpose(syn7)
        xx = [self.repetition7_ste(x) for x in xx]
        return self.transpose(xx)
            
    def syndrome7(self, mat):
        check_rows = [xor(mat[0], mat[i]) for i in range(1,7)]
        return check_rows

    def reset(self):
        self.info.clear()

        
        # check_check = [self.syndrome(x) for x in check_rows] # 6 x 63 (not really)
        # check_rows = [x[63:] for x in check_rows] # 6 x 64
        # check_cols = self.syndrome(mat[0]) # 1 x 63
        
        # return {
        #     "check": check_check,
        #     "rows": check_rows,
        #     "cols": check_cols,
        # }

    def rep_bch_ste(self, syn, err7):
        err = self.syn2error(syn)
        return [xor(err, x) for x in err7]
    
    # def scheme_3(self, r_):
    #     # accumulates seven words in an array
    #     self.info.append(r_[:127])
    #     if len(self.info) < 7:
    #         return [] 
    #     word = list(self.info) 
    #     self.info.clear()

    #     # random code word generator
    #     m = [random.randrange(2) for i in range(64)]
    #     c = self.encode(m)
    #     c = [c for i in range(7)]

    #     # mask with above code word
    #     inp = self.xor2d(c, word)

    #     # calculate syndrome (syn7 + syn)
    #     syn7 = self.syndrome7(inp) # notice that the repetition syndrome that was enrolled will always be 0
    #     syn7 = [self.syndrome(x) for x in syn7]
    #     # err7 = self.stacked_r7ste(syn7)
    #     # syn = self.syndrome(xor(word[0], err7[0]))

    #     # syndrome to error
    #     syn = xor(syn, self.helper["syn7"])
    #     err = self.rep_bch_ste(syn, err7)

    #     # fix errors from original code word
    #     word = self.xor2d(word, err)
    #     #print(word)
    #     return word[0][63:]
    
    
    
    """
    def scheme_3(self, r_):
        # accumulates seven words in an array
        self.info.append(r_[:127])
        if len(self.info) < 7:
            return [] 
        word = list(self.info) 
        self.info.clear()

        # random code word generator
        m = [random.randrange(2) for i in range(64)]
        c = self.encode(m)
        c = [c for i in range(7)]

        # mask with above code word
        inp = self.xor2d(c, word)

        # calculate syndrome (syn7 + syn)
        syn7 = self.syndrome7(inp) # notice that the repetition syndrome that was enrolled will always be 0
        syn = self.syndrome(xor(word[0], err7[0]))

        # syndrome to error
        syn = xor(syn, self.helper["syn7"])
        
        err7 = self.stacked_r7ste(syn7)
        err = self.rep_bch_ste(syn, err7)

        # fix errors from original code word
        word = self.xor2d(word, err)
        #print(word)
        return word[0][63:]
    """

        

    def posture(self, cmd, val, rd):
        val = get_bits(val)
        r = self.bosture(cmd, val)
        return get_bytes(r)
    
    def bosture(self, cmd, val, *vals):
        if self.scheme == None:
            raise Exception("no scheme")
        if cmd == 'e':
            return self.encode(val[:64])
        if cmd == 'p':
            val = val[:127]
            if self.scheme == 1:
                return self.scheme_1(val)
            if self.scheme == 2:
                return self.scheme_2(val)
            if self.scheme == 3:
                return self.scheme_3(val)
        elif cmd == 'k':
            if self.scheme == 1:
                tmp = self.helper["enc"] = xor(self.encode(val[:64]), vals[0])
                return tmp
            if self.scheme == 2:
                tmp = self.helper["syn"] = self.syndrome(val[:127])
                return tmp
            if self.scheme == 3:
                self.info.append(val[:127])
                if len(self.info) < 7:
                    return [] 
                word = list(self.info) 
                self.info.clear()
                syn7 = self.syndrome7(word)
                syn7 = [self.syndrome(x) for x in syn7]
                tmp = self.helper["syn7"] = syn7
                return tmp


    # def bit_trans(self, bits):
    #     bts = []
    #     for i in range(8):
    #        bts.append(bits[i*8:i*8+8])
    #     return bits[56:] + bits[:56]
    
    def encode(self, bits):
        w = get_bytes(bits[::-1])[:8]
        d = self.b.encode(w)
        return get_bits(d)[:63][::-1] + bits

    def decode(self, bits):
        bits = list(bits)
        r = bits[:63][::-1]
        i = bits[63:][::-1]
        r = get_bytes(r)+get_bytes([0])
        i = get_bytes(i)[:8]
        #print(i, r)
        t = self.b.decode(data=i, recv_ecc=r)
        if t<0:
            raise Exception("decoding failed")
        i = get_bits(i)[:64]
        for err in self.b.errloc:
            if err < 64:
                i[err] ^= 1
        return i[::-1]
    
    def syndrome(self, bits):
        bits = list(bits)
        r = bits[:63][::-1]
        i = bits[63:][::-1]
        r = get_bytes(r)+get_bytes([0])
        i = get_bytes(i)[:8]
        # self.b.__dict__["syn"] = 
        t = self.meta["errors"] = self.b.decode(data=i, recv_ecc=r)
        if t==0:
            return tuple(0 for _ in range(20))
        else:
            return self.b.syn#[::-1]
    # the probablem could originate from the fact that only 2*t syndromes are calculated
    def syn2error(self, syn):
        t = self.b.decode(syn=syn)
        error = [0 for _ in range(127)]
        i = [0 for _ in range(64)]
        r = [0 for _ in range(63)]
        for err in self.b.errloc:
            error[err] ^= 1
            continue
            if err>=64:
                i[err-64] ^= 1#error[err-64] ^= 1
            else:
                r[err+63] ^= 1
        return error[::-1]

    def voters(self, syn7):
        xx = self.transpose(syn7)
        xx = [int(sum(x) > 3) for x in xx]
        return xx
    
    def scheme_3(self, r_):
        self.info.append(r_[:127])
        if len(self.info) < 7:
            return [] 

        # print("INSIDE")
        word = list(self.info) 
        self.info.clear()

        # random code word generator
        m = [random.randrange(1) for i in range(64)]
        c = self.encode(m)
        c = [c for i in range(7)]

        # mask with above code word
        inp = self.xor2d(c, word)

        # calculate syndrome (syn7 + syn)
        
        syn7 = self.syndrome7(inp)
        syn7 = [self.syndrome(x) for x in syn7]
        syn7 = self.xor2d(syn7, self.helper["syn7"])
        
        for x in syn7:
            print(x)
            
        err7 = [self.syn2error(x) for x in syn7]
        for x in err7:
            print(x)
        err7 = self.voters(err7)

        print(err7)
        
        word = xor(word[0], err7)
        # print("OUTSIDE")
        return word[63:]



# stmp = Scheme(3)
# ptmp = [[random.randrange(2) for _ in range(127)] for i in range(7)]
# ktmp = ptmp[0][63:]
# for x in ptmp:
#     a = stmp.bosture("k", x)

# for i in range(1,7):
#     for j in range(9):
#         ptmp[i][random.randrange(127)] ^= 1
# for x in ptmp:
#     a = stmp.bosture("p", x)
# print(a == ktmp)
