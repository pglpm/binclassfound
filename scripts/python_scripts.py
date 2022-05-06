def direct_prob(x, cl, w, p, m, s):
    """
    Gives probability of class conditional on RF-output
    x = RF output (between 0 and 1)
    cl = class (0 or 1)
    w = array of weights
    p = array of probs
    m = array of means for gaussian
    s = array of standard devs for gaussian
    """
    c = 1 - 2**(-10)
    x2 = 0.5 + c*(x-0.5)
    y = np.log(x2/(1-x2))
    pc = p*cl + (1-p)*(1-cl)
    pjoint = np.sum(w * pc * norm.pdf(y, loc=m, scale=s))
    px = np.sum(w * norm.pdf(y, loc=m, scale=s))
    pjoint/px


def inverse_prob(x, cl, w, p, m, s):
    """
    Gives probability of RF-output conditional on class
    x = RF output (between 0 and 1)
    cl = class (0 or 1)
    w = array of weights
    p = array of probs
    m = array of means for gaussian
    s = array of standard devs for gaussian
    """
    c = 1 - 2**(-10)
    x2 = 0.5 + c*(x-0.5)
    y = np.log(x2/(1-x2))
    jac = 4*c/(1 - (c*(1-2*x))**2)
    pc = p*cl + (1-p)*(1-cl)
    pjoint = np.sum(w * pc * norm.pdf(y, loc=m, scale=s) * jac)
    pcl = np.sum(w * pc)
    pjoint/pcl

def CNN_direct_prob(x, cl, w, p, m0, s0, m1, s1):
    """
    Gives probability of class conditional on CNN-output
    x = CNN output: vector of two real numbers
    cl = class (0 or 1)
    w = array of weights
    p = array of probs
    m0 = array of means for gaussian, output0
    s0 = array of standard devs for gaussian, output0
    m1 = array of means for gaussian, output1
    s1 = array of standard devs for gaussian, output1
    """
    pc = p*cl + (1-p)*(1-cl)
    pjoint = np.sum(w * pc * norm.pdf(x[0], loc=m0, scale=s0) * norm.pdf(x[1], loc=m1, scale=s1))
    px = np.sum(w * norm.pdf(x[0], loc=m0, scale=s0) * norm.pdf(x[1], loc=m1, scale=s1))
    pjoint/px



def get_expected_utility(um, m):
    um = np.array(um)
    m = np.array(m)
    np.matmul(um, m)



def choose_class(x, um):
    """
    x = either one number between 0-1, or a (normalized) vector of probabilities
        if it's one number then we are in a binary classification case
        and x is the probability of CLASS 0
    um = utility matrix in the format [[T0,F0],[F1,T1]]
    """
    x = np.array(x)
    if len(x)==1: # it means x is the probability of class 0
        x = np.concatenate((x, 1-x)) # transform into prob. vector
    um = np.array(um)
    utilities = np.matmul(um, x)
    np.argmax(utilities)
