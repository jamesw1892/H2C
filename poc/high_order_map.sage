"""
Implementation of Appendices K.3.1 and K.4.1 of the IETF LWIG's Curve
Representations Draft:
https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23

Given a Finite Field F and an Elliptic Curve E over F with non-zero domain
paramters, IE the curve equation is y^2 = f(x) = x^3 + A*x + B for A, B in
F\{0}. Let E(F) be the group of points on E (including the point at infinity).
Let G be the (sub)group of E(F) of largest prime order.

Appendix K.3 describes how to map an element t in F that is not a square in F
to E(F).

Appendix K.4 uses the mapping in appendix K.3 to map an element t in F that is
not a square in F to G.

This replaces the map_to_curve and clear_cofactor stages of the hash-to-curve
draft since it maps directly to G. However it's domain is only non-squares in F.

It also requires the domain parameters to be non-zero, however note 2 in
appendix K.3.1 states that it is often possible to find an isogenous curve
with non-zero domain parameters over the same field F. Then non-square elements
of F can be mapped this isogenous curve and further mapped to the original
curve using the isogeny.
TODO: Look at isogeny maps in H2C appendix E

The parity function in appendix H is also mentioned which assigns a sign to
elements of F - returns 1 iff the element is "negative" in F (so par(0) = 0).
Pseudocode for the same function is provided in the hash-to-curve draft:
https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-14.html#name-the-sgn0-function
And this is implemented in constant time in the `sgn0` function of `common.sage`
so we use this.
"""

import logging
import sys

try:
    from sagelib.common import sgn0 # TODO: do we need big endian sgn0?
except ImportError as e:
    sys.exit("Error loading preprocessed sage files. Try running `make clean pyfiles`. Full error: " + e)

# logging format - log everything and add 'HighOrderMap' to front of messages
logging.basicConfig(level=logging.INFO, format="HighOrderMap %(levelname)s: %(message)s")

class HighOrderMap:
    maps = dict()

    def __init__(self, name: str, F, A, B, p0x, p0y, non_square):
        """
        Name: Curve name
        F: Finite field the curve is defined over (must be element of GF(q) for some q)

        Any of the following can be integers or already elements of F
        A, B: Curve domain parameters such that curve equation is y^2 = x^3 + A*x + B
        p0x, p0y: coordinates of p0 for this curve
        non_square: non-square element of F
        """

        self.name = name
        self.F = F

        # if they are already elements of F this does nothing
        self.A = F(A)
        self.B = F(B)
        self.p0x = F(p0x)
        self.p0y = F(p0y)
        self.non_square = F(non_square)

        # verify curve domain parameters
        assert A != 0 and B != 0, "This mapping does not work when either curve parameter is 0"

        # create elliptic curve object with sage over F with domain parameters A, B
        self.curve = EllipticCurve(F, [A, B])

        # create elliptic curve point object with sage for P0
        self.p0 = self.curve(self.p0x, self.p0y)

        # calculates RHS of curve equation
        self.f = lambda x: x**3 + A*x + B

        # save to static class attribute so it can be easily retrieved by name
        HighOrderMap.maps[name] = self

    def sqrt(self, x):
        """
        Always pick the square root with 0 parity.

        x can be an int or element of F
        """
        assert self.F(x).is_square(), "Must be a square to find sqrt"
        y = self.F(x).sqrt(extend=False) # sage function sqrt
        if sgn0(y) == 1:
            y = -y
        return y

    def k4(self, t):
        """
        https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23#appendix-K.4
        Return a curve point in a prime order subgroup from the element t in the
        finite field

        t can already be an element of F or an integer
        """

        s = 0 # TODO: Does it matter what this binary digit is?

        F = self.F
        A = self.A
        B = self.B
        t = F(t)

        assert A != 0 and B != 0, "This mapping does not work when either curve parameter is 0"
        assert not t.is_square(), "This mapping cannot map elements that are squares in the finite field"

        if t != -1:
            pt = self.k3(t)
            # we know p0, p0 + P(t), p0 - P(t) are all in the large subgroup
            # whatever t is given it is != 1 and non-square in F

            # get x, y coordinates of the point P(t) returned by k3
            try:
                x, y = pt.xy()
            except ZeroDivisionError as e:
                logging.error("Point at infinity obtained from k3!")
                raise e

            # get x, y coordinates of P0
            try:
                x0, y0 = self.p0.xy()
            except ZeroDivisionError as e:
                logging.error("P0 is point at infinity!")
                raise e

            if sgn0(y0*y) == s:
                return self.p0 + pt
            else:
                return self.p0 - pt

        # t = -1 (and -1 is not a square by assert above)
        else:
            return self.p0

    def k3(self, t):
        """
        https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23#appendix-K.3
        Map an element t in F that is not a square in F to a point of a
        Weierstrass curve if neither domain parameter is zero.

        t can already be an element of F or an integer

        Returns a point on the curve
        """

        F = self.F
        A = self.A
        B = self.B
        t = F(t)

        assert A != 0 and B != 0, "This mapping does not work when either curve parameter is 0"
        assert not t.is_square(), "This mapping cannot map elements that are squares in the finite field"

        if t != -1:

            x = (-B/A)*(1+1/(t+t**2))
            fx = self.f(x)
            if fx.is_square():
                return self.curve(x, self.sqrt(fx))
            else:
                x_prime = t*x
                fx_prime = self.f(x_prime)
                return self.curve(x_prime, self.sqrt(fx_prime))

        # t = -1 (and -1 is not a square by assert above)
        else:
            # can map to point at infinity or P0 (if suitable) see note 1
            # return self.curve(0, 1, 0) # TODO: THINK this is identity
            return self.p0 # TODO: suitable?
