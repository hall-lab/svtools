import numpy as np

def get_p(ls):
    '''
    Convert a value in log-space back to probability-space
    '''
    return np.exp(ls)

def get_ls(p):
    '''
    Convert a value probability-space to log-space
    '''
    if p == 0:
        return float("-inf")
    else:
        return np.log(p)

def ls_multiply(x, y):
    '''
    Multiply two values in log-space
    '''
    if (x == float("-inf")) or (y == float("-inf")):
        return float("-inf")
    else:
        return x + y

def ls_divide(x, y):
    '''
    Divide two values in log-space
    '''
    return x - y

def ls_add(x, y):
    '''
    Add two values in log-space
    '''
    if x == float("-inf"):
        return y
    elif y == float("-inf"):
        return x
    # TODO math.log1p may be a better choice or some numpy equivalent
    elif (x < y):
        return y + np.log(1 + np.exp(x - y))
    else:
        return x + np.log(1 + np.exp(y - x))
