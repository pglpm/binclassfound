def RF_direct_prob(x, cl, w, p, m, s):
    """
    Gives probability of class conditional on RF-output
    x = RF output (between 0 and 1)
    cl = class (0 or 1)
    w = array of weights
    p = array of probs
    m = array of means for gaussian
    s = array of standard devs for gaussian
    """
    pc = p*cl + (1-p)*(1-cl)
    pjoint = np.sum(w * pc * norm.pdf(x, loc=m, scale=s))
    px = np.sum(w * norm.pdf(x, loc=m, scale=s))
    pjoint/px


def RF_inverse_prob(x, cl, q, alpha, mu, sigma):
    """
    Gives probability of RF-output conditional on class
    x = RF output (between 0 and 1)
    cl = class (0 or 1)
    q = array of weights
    alpha = array of probs
    mu = array of means for gaussian
    sigma = array of standard devs for gaussian
    """
    pc = alpha*cl + (1-alpha)*(1-cl)
    pjoint = np.sum(q * pc * norm.pdf(x, loc=m, scale=s))
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
