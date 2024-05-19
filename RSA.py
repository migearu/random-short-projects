 # rsa with some stuff :)
# basic rsa explanation
# select 2 prime numbers, (P, Q) usually large prime numbers
# just as an example, we choose P = 137, Q = 193
# public key = P*Q = 26,441
# now we need a small exponent, that is not a factor of the public key
# let us choose e = 5

# quick note, e must be below f(n), and above 1

# -- SO FAR: P = 137, Q = 193, pubkey = 26,441, e = 5 --

# private key time :)
# define f(n), such that f(n) = (P-1)(Q-1)
# f(n) = 136*192 = 26,112
# now the private key can be calculated for any integer k, which results in an integral private key
# prvkey = (k * f(n) + 1) / e, assume k = 2
# prvkey = ((2 * 26112) + 1) / 5 = 10445

# now we have
# public key = 26,441
# private key = 10,445

# -- outline --
# e and f(n) must be coprime
# e must fall within 1 < e < f(n)
# e =/= p, q

import math
import random

# --- important functions crucial to methods --- #

def fastModExp(base, exp, mod):
    r = 1
    if 1 & exp:
        r = base
    while exp != 0:
        exp //= 2
        base = base**2 % mod
        if exp & 1: r = (r * base) % mod
    return r

# --- prime generation --- #

primesList = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
            71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
            149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
            227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
            307, 311, 313, 317, 331, 337, 347, 349]

def nBitRandom(n):
    return random.randrange(2**(n-1)+1, 2**n-1)

def getLowLevelPrime(n):
    while True:
        candidate = nBitRandom(n)
        for divisor in primesList:
            if candidate % divisor == 0 and divisor ** 2 <= candidate:
                break
        else:
            return candidate

def millerRabinTest(candidate):
    maxDivisionsByTwo = 0
    ec = candidate-1
    while ec % 2 == 0:
        ec //= 2
        maxDivisionsByTwo += 1
    assert(2**maxDivisionsByTwo * ec == candidate-1)
    
    def compositeTrial(roundTester):
        if fastModExp(roundTester, ec, candidate) == 1:
            return False
        for i in range(maxDivisionsByTwo):
            if fastModExp(roundTester, 2**i * ec, candidate) == candidate-1:
                return False
        return True
    
    numTrials = len(primesList)
    for i in range(numTrials):
        roundTester = random.randrange(2, candidate)
        if compositeTrial(roundTester):
            return False
    return True

def getPrimeCandidate(n):
    primeCandidate = 0
    millerRabinPassed = False
    while not millerRabinPassed:
        primeCandidate = getLowLevelPrime(n)
        millerRabinPassed = millerRabinTest(primeCandidate)
        #print(primeCandidate, millerRabinPassed)
    return primeCandidate

# --- this is all for rsa --- #

def egcd(a, b): # extended euclidean algorithm
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        (q, a), b = divmod(b, a), a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return b, x0, y0

def modinv(a, b): # this returns x, such that ax mod b = 1
    g, x, _ = egcd(a, b)
    if g != 1:
        return -1
    return x % b

def keyPairGenerator(p, q):
    n = p*q
    fn = math.lcm((p-1), (q-1))
    #print(fn)
    e = 2
    k = 1
    while e < fn:
        if e == p or e == q:
            continue
        if math.gcd(e, fn) == 1:
            break
        e += 1
    d = modinv(e, fn)
    if d < 0: d += fn
    pubkey = (e, n)
    prvkey = (d, n)
    return (pubkey, prvkey)
    
def rsaEncrypt(msg, publicKey):
    msgbytes = str.encode(msg)
    integerRep = int.from_bytes(msgbytes, byteorder='little', signed=False)
    encryptedRep = fastModExp(integerRep, publicKey[0], publicKey[1])
    #print(integerRep)
    #print(encryptedRep.to_bytes(byteorder='little', signed=False).decode('ISO-8859-1'))
    return '{0:x}'.format(encryptedRep)

def rsaDecrypt(enc, privateKey):
    encryptedRep = int(enc, base=16)
    integerRep = fastModExp(encryptedRep, privateKey[0], privateKey[1])
    msgbytes = integerRep.to_bytes(math.ceil(integerRep.bit_length()/8), byteorder='little', signed=False)
    return msgbytes.decode('utf-8')

keyPair = keyPairGenerator(getPrimeCandidate(64), getPrimeCandidate(64))
print(keyPair)

msg = "smash hit"
print(msg)
# say we want to encrypt the message
encrypted = rsaEncrypt(msg, keyPair[0])
print(encrypted)

# now we decrypt our encrypted message
decrypted = rsaDecrypt(encrypted, keyPair[1])
print(decrypted)
