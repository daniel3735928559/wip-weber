## Weber functions and ray class fields

### 20151209 What is this file?

I've had a hard time keeping track of my notebook scraps and random
relevant python files, so in the next few days I'm going to collect
all my existing notes and code here, as much as possible, hopefully
eventually bringing this file up to the current state of my knowledge.
Thereafter, I'll try to use it as a place to store sanitised versions
of my notes.

New entries will be added to the end of this file.

### 20151209 Recap of problem context and statement:

#### Cryptography connections?

Public key crypto algorithms publish a "public key" which can be used
by anyone for encryption, and maintain a secret, related "private key"
which is required to decrypt messages encrypted with the corresponding
public key.  This necessitates that the public key and private key be
related, and in principle will always be possible to deduce the
private from the public.

Effective crypto algorithms are those for which this deduction is not
actually realistic for an attacker with, say, the budget and resources
of a superpower nation-state.  One may further demand that, to be
declared effective, an algorithm must be able to resist an attacker
with a functioning quantum computer.

Standard number theoretic cryptographic algorithms such as ElGamal
(based on the discrete log problem) and RSA (based on the factoring
problem) have been effective to date, assuming the size of the problem
is large enough.  They are both known to be solved by quantum
computers, however.  They also both have the following slightly
unsatisfying feature: Both rely on picking "random prime numbers (one
for ElGamal, and two for RSA) away from a bad set with known attacks".
For example, factoring is relatively easy if the two primes are very
close together.  There is a growing list of features the set of primes
may have that make them easy to factor.

This is fundamentally because it is not known how to reduce the
problem for a given choice of prime(s) to the problem for any other
choice, so it is in principle possible that the problem is easy to
solve in the average case, where the primes are chosen completely at
random, and it is just the worst case (which may in fact be rare)
which is actually hard.

The consequences are twofold:

1. We do not know what set of n-bit primes we can safely pick.

2. We do not even know that increasing n increases the difficulty of
the problem.  For instance, it could be true (though it would be
extremely surprising to my current knowledge) that primes with an even
number of bits are very easy to detect in a factorisation whereas
primes with an odd number of bits are hard.

#### An intermediate problem: Knapsack problems

One attempt at a cryptosystem based on a hard problem is one based on
the "knapsack problem": "Given a set S = {a_1, ..., a_n} and a target
x, find a subset of S whose elements add up to x."  In one
cryptosystem based on this, the private key is a random r and a set S
whose elements are "superincreasing"--that is, `a_{i+1} > a_1 + ... +
a_i` for all i.  (Solving the knapsack problems for superincreasing
sets is easily accomplished by the greedy algorithm.)  The public key
is a modulus m > a_n and the set `rS = {ra_i mod m : i = 1..n}`.  Then
it is easy to encode a message consisting of bits b_1, ..., b_n as the
number sum_i b_i ra_i mod m.  To decode the message, someone knowing
the private key can just multiply by 1/r mod m and solve the
superincreasing knapsack problem trivially.

#### Lattice problems

This system was broken in several different ways, but one observation
is that this is the same as the following problem: Define
(n+1)-dimensional vectors v_i = (a_i, 0, 0, ..., 2, 0, ..., 0) with a
2 in the (i+1)th place, and v_{n+1} = (x, 1, 1, ..., 1).  Then a
solution to the knapsack problem is the same as a linear combination
of these vectors that is within sqrt(n) of the origin.  That is, it is
a solution to the "shortest vector problem" on a specific lattice.

Classic examples of lattice problems include: Given a basis of a
lattice in R^n, finding the shortest vector (SVP), or finding the
shortest n independent vectors (SIVP) or at least finding a vector (or
basis) within a polynomial factor of the shortest.  These problems
have various reductions between them and are generally believed to be
hard in the worst case, much like factoring.  Like factoring, they are
not hard every case, and the lattices constructed in the above way
from superincreasing knapsacks are examples of those susceptible to
one of the well-known algorithms for approimately solving lattice
problems: LLL.

The LLL algorithm solves the SVP within a factor of 2^{(n-1)/2}.  That
is, if l is the length of the actual shortest vector, then LLL will in
polynomial time return a vector of length at most 2^{(n-1)/2} l.  (An
improvement due to Schnorr returns a vector of length at most
(1+epsilon)^n l for any chosen epsilon > 0.)

So once again, it is the case that some lattices have easier solutions
to these "probably hard in the worst case" problems than others.
However, unlike with factoring, the second issue is easily addressed
for lattice problems: If we can quickly solve the shortest vector
problems for lattices in n+1 dimensions, then we can also solve it in
n dimensions by embedding our n vectors v_1, ..., v_n in n+1
dimensions by adding a 0 in the last component, and then letting
v_{n+1} be (0, 0, 0, ..., 999999999999999), (or, generally, some
number much greater than, say, the sum of all the absolute values of
all components of the v_1, ..., v_n).  So at least we are guaranteed
that jacking up the dimension increases the difficulty of these
problems.

The first problem remained, however.  One way of approaching it would
be to generate a class of explicit lattice problems for which a
solution that worked for a random problem from that class would give
rise to a solution to the general lattice problem in the worst case.
This was first accomplished by
[Ajtai](https://www.cs.cmu.edu/~yiwu/doc/ajtai-worst-avg.pdf), who
described a class of problems where the ability to solve a random
instance of the problem implies the ability to polynomially
approximate SIVP in the worst case.

Specifically, Ajtai says for constants A, B, and given n, to pick an
`m > A n log(n)`, and `q = n^B`.  Then random n-dimensional vectors v_1,
..., v_m with integer coordinates chosen uniformly from [0,q-1]
determine a so-called q-ary lattice (a lattice where if two vectors
are congruent mod q then they are either both in the lattice or both
not).  He proves that an algorithm that can find a v_1, ..., v_m just
from knowing the lattice and q can solve approximate SIVP in the worst
case.

#### Problems with lattices for cryptosystems

There are two issues with using this system for cryptography:

1. The function is not injective.  m is a bit larger than n.  This
means the function might be a reasonable hash function, but it is not
on its own necessarily possible to decrypt.

2. Efficiency: As n increases, the amount of data to store grows like
n^2.

We will not talk much about point 1 today (though it is quite
interesting and makes me kind of nervous), but point 2 bears some
remarks.

Micciancio proved a similar result to Ajtai except restricting to only
cyclic lattices: that is, lattices given by taking a A*log(n)
n-dimensional vectors and cycling their coordinates around for m =
A*n*log(n) vectors in total.  This requires only nlog(n) storage,
rather than n^2, which is much better.

However, cyclic lattices no longer have the property that they are
easy to extend one dimension up.  So once again we do not know
whether, for example, prime dimensions of cyclic lattices give hard
problems whereas dimensions with only small prime factors are easy.

#### Ideal lattices

Ideal lattices are a slight generalisation of cyclic lattices which
retain the efficiency property and the average-to-worst-case reduction
(though now "worst case" being in the sense of "worst among ideal
lattices", rather than "worst case among all lattices").  Ideal
lattices come from rings of integers in number fields (or, as it is
often phrased, from irreducible, monic, degree-n polynomials f in Z[x]).


#### Class fields as interesting sources of lattices

It turns out (Pikert and Rosen 07) that one can recover the guarantee
of increasing security from increasing n if one works in a family of
number fields with constant root discriminant.  For example, towers K
= H_0 < H_1 < H_2 < ... where H_i is Hilbert class field of H_{i-1}.
Golod-Shafarevich implies that there are towers of this kind that grow
forever.

So if we can compute Hilbert class fields, then we can write down
these families of lattices and potentially get provably secure crypto
(as always, on the assumption that the worst case of some modified
version of SVP or SIVP is actually hard).

#### Computing class fields

From the theory of complex multiplication on elliptic curves over Q,
we know that the j(E) generates the Hilbert class field of the CM
field of E.  For example, j((d+\sqrt d)/2) generates the Hilbert
class field of k=Q(sqrt d).

However, Zagier and Yui use the following approach: Cl(Q(sqrt d)) is
generated by ideals corresponding to binary quadratic forms of
discriminant d.  They show that there is an SL_2(Z)-invariant function
f on these forms constructed out of the "Weber" functions f_0, f_1,
f_2 (which are permuted with 48th root of unity factors by SL_2(Z), so
it is conceivable one could understand this action enough to construct
something invariant out of them).  The values of this f, then, are the
Galois conjugates of a generator of the Hilbert class field, and as
such are the roots of a polynomial defining this field.  We can
therefore compute this polynomial, and it turns out that the
coefficients we get are interestingly small--much smaller than those
coming from computing the minimal polynomial of j directly.

For example, here is some Sage code:

```
R = ZZ[['q']]
q = R.gen()

qexp_eta(R, 100)

def f(t, prec=100):
    q = e^(2*pi*I*t)
    ans = q^(-1/48)
    for n in range(1,prec):
        ans = ans * (1+q^(n-(1/2)))
    return ans

def f1(t, prec=100):
    q = e^(2*pi*I*t)
    ans = q^(-1/48)
    for n in range(1,prec):
        ans = ans * (1-q^(n-(1/2)))
    return ans

def f2(t, prec=100):
    q = e^(2*pi*I*t)
    ans = sqrt(2)*q^(1/24)
    for n in range(1,prec):
        ans = ans * (1+q^n)
    return ans

def fq(a, b, c, d, prec=100):
    z = e^(2*pi*I/48)
    tQ = (-b + sqrt(d))/(2*a)
    ed = (-1)^((d-1)/8)
    if(a % 2 == 0 and c % 2 == 0):
        return z^(b*(a-c-a*c^2))*f(tQ, prec=prec)
    elif(a % 2 == 0 and c % 2 != 0):
        return ed*z^(b*(a-c-a*c^2))*f1(tQ, prec=prec)
    elif(a % 2 != 0 and c % 2 == 0):
        return ed*z^(b*(a-c+a^2*c))*f2(tQ, prec=prec)
    print("Problem")

def mp(d, prec=100):
    l = BinaryQF_reduced_representatives(d)
    R = PolynomialRing(CC, 'x')
    ans = R(1)
    for Q in l:
        ans *= x - N(fq(Q[0], Q[1], Q[2], d, prec))
    return expand(ans)

#N(fq(1, 1, (1+31)/4, 31, 100))
#N(fq(1, 1, (1+31)/4, 31, 1000))
#N(fq(1, 1, (1+31)/4, 31, 2000))
[round(a.real_part()) for a in PolynomialRing(CC, 'x')(mp(-479,100)).coefficients()]
```

with output

```
1 - q - q^2 + q^5 + q^7 - q^12 - q^15 + q^22 + q^26 - q^35 - q^40 + q^51 + q^57 - q^70 - q^77 + q^92 + O(q^100)
[-1, 10, 23, 76, 104, 126, 109, 144, 205, 317, 336, 280, 138, 11, -76, -82, -48, 6, 47, 66, 60, 41, 22, 9, 3, 1]
```

#### The question

Trying to understand this, f_2 can be written as omega_2^(1/24), which
is a modular form for much bigger group.  Apparently this group is
called \Phi_0^0(24).  Likewise, \omega_2^{1/48} is a modular function
for \Phi_0^0(48), and the associated modular curve has also genus
zero, and its function field is generated by \omega_2^{1/48}. A
question, then is whether \omega_2^{1/48}(\frac{d+\sqrt d}2) generates
the Hilbert class field.  Note that \Phi_0^0(48) is no longer a
congruence subgroup.

