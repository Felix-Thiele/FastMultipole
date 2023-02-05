
def int_bin(val, depth):
    res = str(bin(val))[2:]
    res = '0' * (depth - len(res)) + res
    return res

def bin_int(val):
    if val=='': return 0
    return int(val, 2)
