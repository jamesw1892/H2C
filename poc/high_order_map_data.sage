"""
Stores details from Table 1 of the IETF LWIG's Curve Representations Draft:
https://datatracker.ietf.org/doc/html/draft-ietf-lwig-curve-representations-23
to create instances of the high-order map for different curves. For each curve,
the following attributes are required:

- q (int): The size of the curve's underlying field, must be a prime power
- a (int): The coefficient of x in the curve's Weierstrass equation, must not be 0
- b (int): The constant coefficient in the curve's Weierstrass equation, must not be 0
- h (int): The cofactor of the curve h = |E(F)|/|G| where E(F) is the group of
    points on the elliptic curve and G is the largest subgroup of E(F) with prime order
- P0x (int): The x-coordinate of the point P0 of the curve
- P0y (int or None): The y-coordinate of the point P0 on the curve.
    If None, will be calculated from x
- delta (int or None): A non-square element of the curve's underlying field.
    If None, any such element will be chosen.

P0 is a point on the curve where none of P0, P0 + P(t), nor P0 - P(t) are in
the smallest subgroup of the curve for any non-square element t of the
underlying field where t != -1 and P is the function specified by appendix
K.3.1 of the Curve Representations Draft and implemented in high_order_map.sage.
Sometimes P0 is g, the base point of the curve.

These values are stored and the high-order map instances are created for the
following curves:

In Table 1 of Curve Representations:
- NIST P-224
- NIST P-256
- NIST P-384
- NIST P-521
- Brainpool P224r1
- Brainpool P256r1
- Brainpool P320r1
- Brainpool P384r1
- Brainpool P512r1
- Wei25519
- Wei25519.2
- Wei25519.-3
- Wei448
- Wei448.1
- Wei448.-3
- secp256k1

Not In Table 1 of Curve Representations:
- NIST P-192
- Brainpool P-160 R1
- Brainpool P-192 R1
- BLS12-381 G1 (1)
- BLS12-381 G2 (1, 2)
- secp192k1 (3)
- secp192r1
- secp224k1 (3)
- secp224r1
- secp256r1
- secp384r1
- secp521r1

(1) Haven't found base point
(2) Has unusual field construction
(3) Have a = 0 so require an isogenous curve, but we have not currently found one

Therefore we have not fully defined or instantiated these curves
"""

import sys

try:
    from sagelib.high_order_map import HighOrderMap, IsoHighOrderMap
    from sagelib.suite_p256 import p256_p, p256_A, p256_B
    from sagelib.suite_p384 import p384_p, p384_A, p384_B
    from sagelib.suite_p521 import p521_p, p521_A, p521_B
    from sagelib.suite_secp256k1 import secp256k1_p, Ap as secp256k1_a_iso, Bp as secp256k1_b_iso, iso_map as secp256k1_iso_map
    from sagelib.suite_bls12381g1 import p as bls12381g1_p, Ap as bls12381g1_a_iso, Bp as bls12381g1_b_iso, h_eff as bls12381g1_h, iso_map as bls12381g1_iso_map
    from sagelib.suite_bls12381g2 import p as bls12381g2_p, Ap as bls12381g2_a_iso, Bp as bls12381g2_b_iso, h_eff as bls12381g2_h, iso_map as bls12381g2_iso_map
except ImportError as e:
    sys.exit("Error loading preprocessed sage files. Try running `make clean pyfiles`. Full error: " + e)

###############################################################################
##################### In Table 1 of Curve Representations #####################
###############################################################################

# NIST P-224 - see FIPS-186-4
p224_q = 26959946667150639794667015087019630673557916260026308143510066298881
p224_a = -3
p224_b = 0xb4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4
p224_h = 1
p224_gx = 0xb70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21
p224_gy = 0xbd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34
p224_delta = 11
p224_map = HighOrderMap(p224_q, p224_a, p224_b, p224_h, p224_gx, p224_gy, p224_delta)

# NIST P-256 - see FIPS-186-4
p256_h = 1
p256_P0x = 0
p256_delta = -1
p256_map = HighOrderMap(p256_p, p256_A, p256_B, p256_h, p256_P0x, None, p256_delta)

# NIST P-384 - see FIPS-186-4
p384_h = 1
p384_P0x = 0
p384_delta = -1
p384_map = HighOrderMap(p384_p, p384_A, p384_B, p384_h, p384_P0x, None, p384_delta)

# NIST P-521 - see FIPS-186-4
p521_h = 1
p521_P0x = 0
p521_delta = -1
p521_map = HighOrderMap(p521_p, p521_A, p521_B, p521_h, p521_P0x, None, p521_delta)

# Brainpool P-224 R1 - see RFC 5639
b224_q = 0xD7C134AA264366862A18302575D1D787B09F075797DA89F57EC8C0FF
b224_a = 0x68A5E62CA9CE6C1C299803A6C1530B514E182AD8B0042A59CAD29F43
b224_b = 0x2580F63CCFE44138870713B1A92369E33E2135D266DBB372386C400B
b224_h = 1
b224_gx = 0x0D9029AD2C7E5CF4340823B2A87DC68C9E4CE3174C1E6EFDEE12C07D
b224_gy = 0x58AA56F772C0726F24C6B89E4ECDAC24354B9E99CAA3F6D3761402CD
b224_delta = -1
b224_map = HighOrderMap(b224_q, b224_a, b224_b, b224_h, b224_gx, b224_gy, b224_delta)

# Brainpool P-256 R1 - see RFC 5639
b256_q = 0xA9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377
b256_a = 0x7D5A0975FC2C3057EEF67530417AFFE7FB8055C126DC5C6CE94A4B44F330B5D9
b256_b = 0x26DC5C6CE94A4B44F330B5D9BBD77CBF958416295CF7E1CE6BCCDC18FF8C07B6
b256_h = 1
b256_gx = 0x8BD2AEB9CB7E57CB2C4B482FFC81B7AFB9DE27E1E3BD23C23A4453BD9ACE3262
b256_gy = 0x547EF835C3DAC4FD97F8461A14611DC9C27745132DED8E545C1D54C72F046997
b256_delta = -1
b256_map = HighOrderMap(b256_q, b256_a, b256_b, b256_h, b256_gx, b256_gy, b256_delta)

# Brainpool P-320 R1 - see RFC 5639
b320_q = 0xD35E472036BC4FB7E13C785ED201E065F98FCFA6F6F40DEF4F92B9EC7893EC28FCD412B1F1B32E27
b320_a = 0x3EE30B568FBAB0F883CCEBD46D3F3BB8A2A73513F5EB79DA66190EB085FFA9F492F375A97D860EB4
b320_b = 0x520883949DFDBC42D3AD198640688A6FE13F41349554B49ACC31DCCD884539816F5EB4AC8FB1F1A6
b320_h = 1
b320_gx = 0x43BD7E9AFB53D8B85289BCC48EE5BFE6F20137D10A087EB6E7871E2A10A599C710AF8D0D39E20611
b320_gy = 0x14FDD05545EC1CC8AB4093247F77275E0743FFED117182EAA9C77877AAAC6AC7D35245D1692E8EE1
b320_delta = -1
b320_map = HighOrderMap(b320_q, b320_a, b320_b, b320_h, b320_gx, b320_gy, b320_delta)

# Brainpool P-384 R1 - see RFC 5639
b384_q = 0x8CB91E82A3386D280F5D6F7E50E641DF152F7109ED5456B412B1DA197FB71123ACD3A729901D1A71874700133107EC53
b384_a = 0x7BC382C63D8C150C3C72080ACE05AFA0C2BEA28E4FB22787139165EFBA91F90F8AA5814A503AD4EB04A8C7DD22CE2826
b384_b = 0x04A8C7DD22CE28268B39B55416F0447C2FB77DE107DCD2A62E880EA53EEB62D57CB4390295DBC9943AB78696FA504C11
b384_h = 1
b384_gx = 0x1D1C64F068CF45FFA2A63A81B7C13F6B8847A3E77EF14FE3DB7FCAFE0CBD10E8E826E03436D646AAEF87B2E247D4AF1E
b384_gy = 0x8ABE1D7520F9C2A45CB1EB8E95CFD55262B70B29FEEC5864E19C054FF99129280E4646217791811142820341263C5315
b384_delta = -1
b384_map = HighOrderMap(b384_q, b384_a, b384_b, b384_h, b384_gx, b384_gy, b384_delta)

# Brainpool P-512 R1 - see RFC 5639
b512_q = 0xAADD9DB8DBE9C48B3FD4E6AE33C9FC07CB308DB3B3C9D20ED6639CCA703308717D4D9B009BC66842AECDA12AE6A380E62881FF2F2D82C68528AA6056583A48F3
b512_a = 0x7830A3318B603B89E2327145AC234CC594CBDD8D3DF91610A83441CAEA9863BC2DED5D5AA8253AA10A2EF1C98B9AC8B57F1117A72BF2C7B9E7C1AC4D77FC94CA
b512_b = 0x3DF91610A83441CAEA9863BC2DED5D5AA8253AA10A2EF1C98B9AC8B57F1117A72BF2C7B9E7C1AC4D77FC94CADC083E67984050B75EBAE5DD2809BD638016F723
b512_h = 1
b512_gx = 3
b512_delta = -1
b512_map = HighOrderMap(b512_q, b512_a, b512_b, b512_h, b512_gx, None, b512_delta)

# Wei25519 - see LWIG Curve Representations Appendices E.3 and G.3
wei25519_q = 2^255 - 19
wei25519_a = 19298681539552699237261830834781317975544997444273427339909597334573241639236
wei25519_b = 55751746669818908907645289078257140818241103727901012315294400837956729358436
wei25519_h = 8
wei25519_P0x = 3
wei25519_delta = 2
wei25519_map = HighOrderMap(wei25519_q, wei25519_a, wei25519_b, wei25519_h, wei25519_P0x, None, wei25519_delta)

wei25519_2_a = 2
wei25519_2_b = 12102640281269758552371076649779977768474709596484288167752775713178787220689
wei25519_2_P0x = 244
wei25519_2_map = HighOrderMap(wei25519_q, wei25519_2_a, wei25519_2_b, wei25519_h, wei25519_2_P0x, None, wei25519_delta)

wei25519_n3_a = -3
wei25519_n3_b = 29689592517550930188872794512874050362622433571298029721775200646451501277098
wei25519_n3_P0x = 41
wei25519_n3_map = HighOrderMap(wei25519_q, wei25519_n3_a, wei25519_n3_b, wei25519_h, wei25519_n3_P0x, None, wei25519_delta)

# Wei448 - see LWIG Curve Representations Appendices M.3 and N.3
wei448_q = 2^448 - 2^224 - 1
wei448_a = 484559149530404593699549205258669689569094240458212040187660132787074885444487181790930922465784363953392589641229091574035657199637535
wei448_b = 269199527516891440944194002921483160871719022476784466770922295992819380802492878772739401369880202196329216467349495319191685664513904
wei448_h = 4
wei448_P0x = 18
wei448_delta = -1
wei448_map = HighOrderMap(wei448_q, wei448_a, wei448_b, wei448_h, wei448_P0x, None, wei448_delta)

wei448_1_a = 1
wei448_1_b = 659612817018071705319448049859079902872252480565600363923809459513818308850763543778602104492771511922449740791489579066934526889652743
wei448_1_P0x = 10
wei448_1_map = HighOrderMap(wei448_q, wei448_1_a, wei448_1_b, wei448_h, wei448_1_P0x, None, wei448_delta)

wei448_n3_a = -3
wei448_n3_b = 699937686810001500848336699619005330673833355924944987095346934649131425073158306877468995089322968102492731574779458733142208859254465
wei448_n3_P0x = 8
wei448_n3_map = HighOrderMap(wei448_q, wei448_n3_a, wei448_n3_b, wei448_h, wei448_n3_P0x, None, wei448_delta)

# secp256k1 - see SEC 2
secp256k1_h = 1
secp256k1_P0x = 0
secp256k1_delta = -1
secp256k1_map = IsoHighOrderMap(secp256k1_iso_map, secp256k1_p, secp256k1_a_iso, secp256k1_b_iso, secp256k1_h, secp256k1_P0x, None, secp256k1_delta)

###############################################################################
################### Not In Table 1 of Curve Representations ###################
###############################################################################

# NIST P-192 - see FIPS-186-4
p192_q = 6277101735386680763835789423207666416083908700390324961279
p192_a = -3
p192_b = 0x64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1
p192_h = 1
p192_P0x = 0
p192_delta = -1
p192_map = HighOrderMap(p192_q, p192_a, p192_b, p192_h, p192_P0x, None, p192_delta)

# Brainpool P-160 R1 - see RFC 5639
b160_q = 0xE95E4A5F737059DC60DFC7AD95B3D8139515620F
b160_a = 0x340E7BE2A280EB74E2BE61BADA745D97E8F7C300
b160_b = 0x1E589A8595423412134FAA2DBDEC95C8D8675E58
b160_h = 1
b160_gx = 0xBED5AF16EA3F6A4F62938C4631EB5AF7BDBCDBC3 # TODO: is this suitable for P0?
b160_gy = 0x1667CB477A1A8EC338F94741669C976316DA6321 # TODO: is this suitable for P0?
b160_delta = -1
b160_map = HighOrderMap(b160_q, b160_a, b160_b, b160_h, b160_gx, b160_gy, b160_delta)

# Brainpool P-192 R1 - see RFC 5639
b192_q = 0xC302F41D932A36CDA7A3463093D18DB78FCE476DE1A86297
b192_a = 0x6A91174076B1E0E19C39C031FE8685C1CAE040E5C69A28EF
b192_b = 0x469A28EF7C28CCA3DC721D044F4496BCCA7EF4146FBF25C9
b192_h = 1
b192_gx = 0xC0A0647EAAB6A48753B033C56CB0F0900A2F5C4853375FD6 # TODO: is this suitable for P0?
b192_gy = 0x14B690866ABD5BB88B5F4828C1490002E6773FA2FA299B8F # TODO: is this suitable for P0?
b192_delta = -1
b192_map = HighOrderMap(b192_q, b192_a, b192_b, b192_h, b192_gx, b192_gy, b192_delta)

# BLS12-381 G1 - see https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-pairing-friendly-curves-10
# bls12381g1_P0x = 0 # TODO: find
# bls12381g1_P0y = 0 # TODO: find
bls12381g1_delta = -1
# bls12381g1_map = IsoHighOrderMap(bls12381g1_iso_map, bls12381g1_p, bls12381g1_a_iso, bls12381g1_b_iso, bls12381g1_h, bls12381g1_P0x, bls12381g1_P0y, bls12381g1_delta)

# BLS12-381 G2 - see https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-pairing-friendly-curves-10
# bls12381g2_P0x = 0 # TODO: find
# bls12381g2_P0y = 0 # TODO: find
# bls12381g2_delta = None # TODO: find
# bls12381g2_map = IsoHighOrderMap(bls12381g2_iso_map, bls12381g2_p, bls12381g2_a_iso, bls12381g2_b_iso, bls12381g2_h, bls12381g2_P0x, bls12381g2_P0y, bls12381g2_delta)

# secp192k1 - see SEC 2 - TODO: Find isogenous curve
# secp192k1_iso_map
secp192k1_q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFEE37
# secp192k1_a_iso
# secp192k1_b_iso
secp192k1_h = 1
secp192k1_gx = 0xDB4FF10EC057E9AE26B07D0280B7F4341DA5D1B1EAE06C7D # TODO: is this suitable for P0?
secp192k1_gy = 0x9B2F2F6D9C5628A7844163D015BE86344082AA88D95E2F9D # TODO: is this suitable for P0?
secp192k1_delta = -1
# secp192k1_map = IsoHighOrderMap(secp192k1_iso_map, secp192k1_q, secp192k1_a_iso, secp192k1_b_iso, secp192k1_h, secp192k1_gx, secp192k1_gy, secp192k1_delta)

# secp192r1 - see SEC 2
secp192r1_q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF
secp192r1_a = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC
secp192r1_b = 0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1
secp192r1_h = 1
secp192r1_P0x = 0
secp192r1_delta = -1
secp192r1_map = HighOrderMap(secp192r1_q, secp192r1_a, secp192r1_b, secp192r1_h, secp192r1_P0x, None, secp192r1_delta)

# secp224k1 - see SEC 2 - TODO: Find isogenous curve
# secp224k1_iso_map
secp224k1_q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFE56D
# secp224k1_a_iso
# secp224k1_b_iso
secp224k1_h = 1
secp224k1_gx = 0xA1455B334DF099DF30FC28A169A467E9E47075A90F7E650EB6B7A45C # TODO: is this suitable for P0?
secp224k1_gy = 0x7E089FED7FBA344282CAFBD6F7E319F7C0B0BD59E2CA4BDB556D61A5 # TODO: is this suitable for P0?
secp224k1_delta = 2
# secp224k1_map = IsoHighOrderMap(secp224k1_iso_map, secp224k1_q, secp224k1_a_iso, secp224k1_b_iso, secp224k1_h, secp224k1_gx, secp224k1_gy, secp224k1_delta)

# secp224r1 - see SEC 2
secp224r1_q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001
secp224r1_a = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE
secp224r1_b = 0xB4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4
secp224r1_h = 1
secp224r1_gx = 0xB70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21 # TODO: is this suitable for P0?
secp224r1_gy = 0xBD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34 # TODO: is this suitable for P0?
secp224r1_delta = 11
secp224r1_map = HighOrderMap(secp224r1_q, secp224r1_a, secp224r1_b, secp224r1_h, secp224r1_gx, secp224r1_gy, secp224r1_delta)

# secp256r1 - see SEC 2
secp256r1_q = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
secp256r1_a = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC
secp256r1_b = 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B
secp256r1_h = 1
secp256r1_gx = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296 # TODO: is this suitable for P0?
secp256r1_gy = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5 # TODO: is this suitable for P0?
secp256r1_delta = -1
secp256r1_map = HighOrderMap(secp256r1_q, secp256r1_a, secp256r1_b, secp256r1_h, secp256r1_gx, secp256r1_gy, secp256r1_delta)

# secp384r1 - see SEC 2
secp384r1_q = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF
secp384r1_a = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFC
secp384r1_b = 0xB3312FA7E23EE7E4988E056BE3F82D19181D9C6EFE8141120314088F5013875AC656398D8A2ED19D2A85C8EDD3EC2AEF
secp384r1_h = 1
secp384r1_P0x = 0
secp384r1_delta = -1
secp384r1_map = HighOrderMap(secp384r1_q, secp384r1_a, secp384r1_b, secp384r1_h, secp384r1_P0x, None, secp384r1_delta)

# secp521r1 - see SEC 2
secp521r1_q = 0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
secp521r1_a = 0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC
secp521r1_b = 0x0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00
secp521r1_h = 1
secp521r1_gx = 0x00C6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66 # TODO: is this suitable for P0?
secp521r1_gy = 0x011839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650 # TODO: is this suitable for P0?
secp521r1_delta = -1
secp521r1_map = HighOrderMap(secp521r1_q, secp521r1_a, secp521r1_b, secp521r1_h, secp521r1_gx, secp521r1_gy, secp521r1_delta)
