# High-Order Construction Implementation

This extension of the CFRG reference Hash to Curve proof-of-concept implements the High Order Constructions presented in the Curve Representations LWIG draft.

Aside from this README, the only files we have edited are `high_order_map.sage` and `high_order_map_data.sage`. Their usage is explained below.

# Links

- [Original Reference Hash to Curve Proof-of-Concept](https://github.com/cfrg/draft-irtf-cfrg-hash-to-curve/tree/main/poc)
- [This extended Proof-of-Concept](https://github.com/jamesw1892/H2C/tree/main/poc)
- [Hash to Curve CFRG Draft Specification](https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/)
- [Curve Representations LWIG Draft Specification](https://datatracker.ietf.org/doc/draft-ietf-lwig-curve-representations/)

# Usage

Before first use, test the existing proof-of-concept by the CFRG including generating test vectors to make sure everything is working correctly with `make test vectors`.

Our additions to the code do not do anything as-is, instead they are intended to be used in other code:

- The `HighOrderMap` and `IsoHighOrderMap` classes defined in `high_order_map.sage` can be instantiated with the necessary attributes to define the high order map. Then any of the high-order constructions can be performed on any input by running any of the following methods:
    - `k3`
    - `k4`
    - `k5`
    - `k6`
    - `uniform`
    - Also `pickDelta` and `verifyDelta` can be used to calculate `delta` for the curve and check a particular delta is valid respectively
- The constants defined in `high_order_map_data.sage` can be imported to be used elsewhere. If `name` is a curve name, then the following constants are provided (except where specified).
    - `name_q`: The order of the underlying field
    - `name_a`: The coefficient of x in the elliptic curve equation
    - `name_b`: The constant coefficient in the elliptic curve equation
    - `name_h`: The cofactor of the curve
    - `name_P0x`: The x-coordinate of the point P0 required by the high-order constructions. In cases where this is equal to the base point of the curve, this is called `name_gx`. If known for that curve
    - `name_P0y`: The y-coordinate of the point P0 required by the high-order constructions. In cases where this is equal to the base point of the curve, this is called `name_gy`. This can be efficiently calculated from the corresponding x-coordinate so is often not provided. If known for that curve
    - `name_delta`: A non-square element of the underlying field - if known for that curve
- Using the instances of `HighOrderMap` and `IsoHighOrderMap` that are created in `high_order_map_data.sage` in the same way as explained in the first bullet point. Note that some curves that we have found delta but not P0 for are instantiated with their base point despite us not knowing whether this is a suitable P0. The code states clearly when this is the case