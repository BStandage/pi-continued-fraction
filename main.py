from decimal import Decimal as Dec
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns


# experimental_guess:d
# accept: a list of coefficients and its length
# return: the probability distribution of coefficients
def experimental_guess(coef):
    distribution = np.zeros(max(coefficients))
    n = len(coefficients)


    # enumerate through all coefficients of the continued fraction
    # increment distribution at each index corresponding to the coefficients seen
    for i in range(0, n):
        distribution[coef[i] - 1] += 1

    for i in range(len(distribution)):
        distribution[i] = Dec(distribution[i])/Dec(n)

    return distribution


def continued_fraction(n, d):
        dict = {}
        coefficients = []
        while d:
            q = n//d
            coefficients.append(q)
            if q in dict:
                dict[q] += 1
            else:
                dict[q] = 1
            n, d = d, n % d
        print(coefficients)
        return (coefficients, dict)




def sqrt(n, one):
    """
    Return the square root of n as a fixed point number with the one
    passed in.  It uses a second order Newton-Raphson convgence.  This
    doubles the number of significant figures on each iteration.
    """
    # Use floating point arithmetic to make an initial guess
    floating_point_precision = 10**16
    n_float = float((n * floating_point_precision) // one) / floating_point_precision
    x = (int(floating_point_precision * math.sqrt(n_float)) * one) // floating_point_precision
    n_one = n * one
    while 1:
        x_old = x
        x = (x + n_one // x) // 2
        if x == x_old:
            break
    return x

def pi_chudnovsky_bs(digits):
    """
    Compute int(pi * 10**digits)

    This is done using Chudnovsky's series with binary splitting
    """
    C = 640320
    C3_OVER_24 = C**3 // 24
    def bs(a, b):
        """
        Computes the terms for binary splitting the Chudnovsky infinite series

        a(a) = +/- (13591409 + 545140134*a)
        p(a) = (6*a-5)*(2*a-1)*(6*a-1)
        b(a) = 1
        q(a) = a*a*a*C3_OVER_24

        returns P(a,b), Q(a,b) and T(a,b)
        """
        if b - a == 1:
            # Directly compute P(a,a+1), Q(a,a+1) and T(a,a+1)
            if a == 0:
                Pab = Qab = 1
            else:
                Pab = (6*a-5)*(2*a-1)*(6*a-1)
                Qab = a*a*a*C3_OVER_24
            Tab = Pab * (13591409 + 545140134*a) # a(a) * p(a)
            if a & 1:
                Tab = -Tab
        else:
            # Recursively compute P(a,b), Q(a,b) and T(a,b)
            # m is the midpoint of a and b
            m = (a + b) // 2
            # Recursively calculate P(a,m), Q(a,m) and T(a,m)
            Pam, Qam, Tam = bs(a, m)
            # Recursively calculate P(m,b), Q(m,b) and T(m,b)
            Pmb, Qmb, Tmb = bs(m, b)
            # Now combine
            Pab = Pam * Pmb
            Qab = Qam * Qmb
            Tab = Qmb * Tam + Pam * Tmb
        return Pab, Qab, Tab
    # how many terms to compute
    DIGITS_PER_TERM = math.log10(C3_OVER_24/6/2/6)
    N = int(digits/DIGITS_PER_TERM + 1)
    # Calclate P(0,N) and Q(0,N)
    P, Q, T = bs(0, N)
    one = 10**digits
    sqrtC = sqrt(10005*one, one)
    return (Q*426880*sqrtC) // T


def check_accuracy(test, guess):
    for i, (t, s) in enumerate(zip(test, guess)):
        if t != s:
            return i
    return None


def lcm(a, b):
    return abs(a*b) // math.gcd(a, b)


def difference(a, b):
    d = lcm(a[1], b[1])
    numa = d / a[1] * a[0]
    numb = d / b[1] * b[0]
    n = numa - numb

    return n, d



def approximate_Khinchin(coeffs):
    product = 1
    for c in coeffs:
        product *= c

    return product**(1/len(coeffs))


if __name__ == '__main__':
    # select number of decimal places you wish to generate
    n = 10000
    pi = pi_chudnovsky_bs(n)
    accurate_pi = pi_chudnovsky_bs(100000)
    print('Pi to 100,000 decimal places: ' + str(accurate_pi))
    print('Pi to 10,000 decimal places: ' + str(pi))

    #print(int(str(accurate_pi)[10000:]))
    diff = '.' + '0'*10000 + str(accurate_pi)[10000:]
    print('Difference: ' + diff)

    # read txt file of coefficients of pi
    f = open('pi_coefficients.txt', 'r')
    test_coeffs = []
    [test_coeffs.append(x.split(' ')[1]) for x in f.readlines()]
    f.close()



    # convert to rational
    pi = str(pi)[0] + '.' + str(pi)[1:]
    rational = Dec(pi).as_integer_ratio()

    accurate_pi = str(accurate_pi)[0] + '.' + str(accurate_pi)[1:]
    r_a_pi = Dec(accurate_pi).as_integer_ratio()
    #diff = difference(r_a_pi, rational)
    #print(diff[0])
    #print(diff[1])

    # compute coefficeints of continued fraction
    results = continued_fraction(rational[0], rational[1])
    coefficients = results[0]
    frequency_dict = results[1]
    dist = experimental_guess(coefficients)

    print('Wth a decimal length of ' + str(n) + ' there are a total of ' + str(len(coefficients)) + ' coefficients.')
    accuracy = check_accuracy(test_coeffs, coefficients)
    if accuracy == 0:
        print('All ' + str(len(coefficients)) + ' coefficients are accurate.')
    else:
        print('Only ' + str(len(coefficients) - int(accuracy)) + ' coefficients computed correctly.')


    # plot a histogram of the coefficients
    #plt.hist(coefficients)
    #sns.distplot(coefficients)

    expectations = []
    for i in range(1, len(coefficients)):
        expectations.append((1/math.log(2)) * (2 * math.log(1 + i) - math.log(i) - math.log(2+i)))
    print(expectations[:10])
    print(dist[:10])




    sns.lineplot(range(1,len(expectations)+1), expectations, color='r', label='Expectation')
    sns.lineplot(range(1, len(dist) + 1), dist, color='b', label='Guess')
    plt.legend(prop={'size': 11})
    plt.xlim(0, 10)
    plt.show()
    print('Approximation of Khinchin\'s constant: ' + str(approximate_Khinchin(coefficients[1:725])))

    ## TO DO ##

    # show that approximation of pi is accurate
    # show the speed of calculating generic vs non generic continued fractions vs some common functions (e.g. log n, n, n^k, e^x, etc...)
    # draw graph of all speeds
    # draw graph of number of decimals vs accurately generated coefficeints
    # overlay histograms of frequencies of generic vs non generic coefficients (or density plots)
    # compare to expected distribution of pi given on 3/10
