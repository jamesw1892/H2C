"""
Implementation of Appendices K.3.1, K.4.1, K.5 and K.6 of the IETF LWIG's Curve
Representations Draft:
https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23

Given a integer field size q, and non-zero integers a and b, we define the
finite field F of size q and the elliptic curve E over F with Weierstrass
equation y^2 = f(x) = x^3 + a*x + b where a and b are transformed into elements
of F. Let E(F) be the group of points on E (including the point at infinity).

Appendix K.3 describes how to map an element t in F where t is not a square in
F to some point in E(F) (not necessarily a high-order point).

Appendix K.4 describes how to use K.3 to map an element t in F where t is not a
square in F to some high-order point in E(F).

Appendix K.6 describes how to extend the mappings in K.3 or K.4 to take as
input any element u in F (not necessarily a non-square).

Only about 3/8 of the points in E(F) can be mapped to with the above mappings.
Appendix K.5 describes how to extend them to produce a distribution that is
statistically indistinguishable from the uniform distribution over E(F).
However, this loses guarantee of K.4 to not map to a point in a small subgroup.

These require the domain parameters to be non-zero, however note 2 in appendix
K.3.1 states that it is often possible to find an isogenous curve with non-zero
domain parameters over the same field F. Then non-square elements of F can be
mapped this isogenous curve and further mapped to the original curve using the
isogeny.
TODO: Look at isogeny maps in H2C appendix E

The parity function in appendix H is also mentioned which assigns a sign to
elements of F - returns 1 iff the element is "negative" in F (so par(0) = 0).
Pseudocode for the same function is provided in the hash-to-curve draft:
https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-14.html#name-the-sgn0-function
And this is implemented in constant time in the `sgn0` function of
`common.sage` so we use this.
"""

import logging
import sys

try:
    from sagelib.common import sgn0 # TODO: do we need big endian sgn0?
except ImportError as e:
    sys.exit("Error loading preprocessed sage files. Try running `make clean pyfiles`. Full error: " + e)

# logging format - log everything except debug and add 'HighOrderMap' to front of messages
logging.basicConfig(level=logging.INFO, format="HighOrderMap %(levelname)s: %(message)s")

class HighOrderMap:
    """
    GENERAL ATTRIBUTES:

    q (int): Field size, prime power
    a (int): Coefficient of x in Weierstrass curve equation
    b (int): Constant in Weierstrass curve equation
    h (int): Cofactor of the curve
    F: Finite Field object from Sage created with GF(q)
    curve: Elliptic Curve object from Sage created with EllipticCurve(.)
    f (int -> int): Calculates the RHS of the Weierstrass curve equation given x
    identity: Elliptic curve point object from Sage that is the point at infinity,
    the identity element of the group of points

    HIGH ORDER CONSTRUCTION ATTRIBUTES:

    delta (F): A non-square element of F
    P0: A point in curve where none of P0, P0+P(t), nor P0-P(t) are in the
    smallest subgroup for any non-square t != -1 and where P is K.3.1
    P0x (F): x coordinate of P0, element of the field F
    P0y (F): y coordinate of P0, element of the field F
    P1: A point in the largest subgroup of the curve
    P2: A point in the largest subgroup of the curve
    """

    def __init__(self, q: int, a: int, b: int, h: int, P0x: int, P0y: int=None, delta: int=None):
        """
        Calculate things about the curve with field size q, Weierstrass
        curve coefficients a and b, and cofactor h.

        q must be a prime power otherwise ValueError will be thrown.

        NOTE: the high order constructions require a != 0 and b != 0 so instead
        use an isogenous curve and input its a and b here, then map the output
        point to the original curve using the isogeny.

        - P0x: the x coordinate of P0

        Optionally provide (if not provided, will calculate):
        - delta: a non-square element of the of size q
        - P0y: the y coordinate of P0

        Where neither P0, P0 + P(t), nor P0 - P(t) are in a small subgroup for
        any non-square t in the field of size q and P(t) is the K.3 mapping.
        """

        self.q = q
        self.F = GF(q)
        self.a = self.F(a)
        self.b = self.F(b)

        # verify curve domain parameters
        assert self.a != 0 and self.b != 0, "This mapping does not work when either curve parameter is 0"

        self.h = h
        self.curve = EllipticCurve(self.F, [a, b])

        self.identity = self.curve(0, 1, 0)

        # curve equation RHS
        self.f = lambda x: x**3 + self.a * x + self.b

        # calculate delta if not provided
        if delta is None:
            self.delta = self.pickDelta()

        # otherwise check valid
        else:
            self.delta = self.F(delta)
            self.verifyDelta(self.delta)

        self.P0x = self.F(P0x)

        # if y coordinate not provided, calculate from x
        if P0y is None:
            P0y = self.sqrt(self.f(P0x))

        self.P0y = self.F(P0y)
        self.P0 = self.curve(self.P0x, self.P0y)

        # for now set all equal, TODO: potentially change
        self.P2 = self.P1 = self.P0

    def pickDelta(self):
        """
        Return a non-square element of the field
        """

        # about half are non-square so certainly exists
        # brute force
        for e in self.F:
            if not e.is_square():
                return e

    def sqrt(self, x):
        """
        Always pick the even square root.

        x can be an int or element of F
        """

        # make x an element of F
        x = self.F(x)

        assert x.is_square(), "Must be a square to find sqrt"
        y = x.sqrt(extend=False) # sage function sqrt

        # ensure we pick the even sqrt
        if sgn0(y) == 1:
            y = -y

        return y

    def k3(self, t):
        """
        https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23#appendix-K.3
        Return a curve point from an element t in F that is not a square in F.
        Neither domain parameter in the curve's Weierstrass equation can be zero.

        t can already be an element of F or an integer
        """

        t = self.F(t)

        assert not t.is_square(), "This mapping cannot map elements that are squares in the finite field"

        # K.4 never calls with t = -1
        if t == -1:
            return self.identity

        # unique solution to f(t*x) = t^3 * f(x)
        x = -(self.b/self.a)*(1+1/(t+t**2))

        fx = self.f(x)
        if fx.is_square():
            y = self.sqrt(fx)

        # if f(x) is not a square then f(t*x) is
        else:
            x = t*x
            y = -self.sqrt(self.f(x))

        # create curve point object
        return self.curve(x, y)

    def k4(self, t, s):
        """
        https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23#appendix-K.4
        Return a high-order curve point from the non-square element
        t in F and binary digit s.

        Note this could be in the largest subgroup, or another coset. It has
        either the largest prime order or the same order as the whole group.

        t can already be an element of F or an integer
        """

        assert s == 0 or s == 1, "s must be a binary digit"

        t = self.F(t)

        if t == -1:
            return self.P0

        Pt = self.k3(t)
        # we know P0, P0 + P(t), P0 - P(t) are not in a small subgroup
        # whatever t is given it is != 1 and non-square in F

        # get x, y coordinates of the point P(t) returned by k3
        try:
            Ptx, Pty = Pt.xy()
        except ZeroDivisionError as e:
            logging.error("Point at infinity obtained from k3!")
            raise e

        if sgn0(self.P0y * Pty) == s:
            return self.P0 + Pt
        else:
            return self.P0 - Pt

    def k6(self, u, s):
        """
        https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23#appendix-K.6
        Return a curve point not in the small subgroup from the element u in F
        and binary digit s.

        Note this could be in the largest subgroup, or another coset. It has
        either the largest prime order or the same order as the whole group.

        u can already be an element of F or an integer.

        This uses K.4 (not just K.3)
        """

        u = self.F(u)

        if u == 0:
            return self.P1

        return self.k4(self.delta * u**2, s)

    def k5(self, u1, s1, u2, s2):
        """
        https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23#appendix-K.5
        Return a uniform curve point from elements u1, u2 in F and binary digits
        s1, s2.
        
        Note that despite using K.6 with K.4, this can produce a point in the
        small subgroup because the points produced by k6 may be in a coset of
        the largest prime subgroup rather than it itself.

        u1 and u2 can already be an element of F or integers.

        This uses K.6 (the version that uses K.4)
        """

        return self.k6(u1, s1) - self.k6(u2, s2)

    def verifyDelta(self, delta):
        assert not delta.is_square(), "Delta is a square in F!"
