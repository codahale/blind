use curve25519_dalek::constants::{RISTRETTO_BASEPOINT_POINT, RISTRETTO_BASEPOINT_TABLE};
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::MultiscalarMul;
use digest::consts::U64;
use digest::Digest;
use rand::thread_rng;

pub fn key_gen() -> ((Scalar,), (RistrettoPoint, RistrettoPoint)) {
    let x = Scalar::random(&mut rand::thread_rng());
    let x_p = &RISTRETTO_BASEPOINT_TABLE * &x;
    let z = RistrettoPoint::random(&mut thread_rng());

    ((x,), (x_p, z))
}

pub fn sign_s1(
    sk: (Scalar,),
    pk: (RistrettoPoint, RistrettoPoint),
) -> ((Scalar, Scalar, Scalar, Scalar), (RistrettoPoint, RistrettoPoint)) {
    let (x,) = sk;
    let (_, z) = pk;

    let a = Scalar::random(&mut rand::thread_rng());
    let t = Scalar::random(&mut rand::thread_rng());
    let y = Scalar::random(&mut rand::thread_rng());

    let a_p = &RISTRETTO_BASEPOINT_TABLE * &a;
    let c = RistrettoPoint::multiscalar_mul([&t, &y], [&G, &z]);

    ((a, y, t, x), (a_p, c))
}

pub fn sign_s2(sts: (Scalar, Scalar, Scalar, Scalar), c: Scalar) -> (Scalar, Scalar, Scalar) {
    let (a, y, t, x) = sts;
    let s = (c * y * x) + a;
    (s, y, t)
}

pub fn sign_u1<D>(
    pk: (RistrettoPoint, RistrettoPoint),
    msg1: (RistrettoPoint, RistrettoPoint),
    m: &[u8],
) -> (
    (
        Scalar,
        Scalar,
        Scalar,
        Scalar,
        Scalar,
        Scalar,
        RistrettoPoint,
        RistrettoPoint,
        RistrettoPoint,
        Scalar,
    ),
    Scalar,
)
where
    D: Digest<OutputSize = U64>,
{
    let (x, z) = pk;
    let (a, c) = msg1;

    let r1 = Scalar::random(&mut rand::thread_rng());
    let r2 = Scalar::random(&mut rand::thread_rng());
    let gamma1 = Scalar::random(&mut rand::thread_rng());
    let gamma2 = Scalar::random(&mut rand::thread_rng());

    let a_p = RistrettoPoint::multiscalar_mul([&r1, &(gamma1 * gamma2.invert())], [&G, &a]);
    let c_p = RistrettoPoint::multiscalar_mul([&gamma1, &r2], [&c, &G]);

    let d = D::new()
        .chain_update(a_p.compress().as_bytes())
        .chain_update(c_p.compress().as_bytes())
        .chain_update(m);

    let mut output = [0u8; 64];
    output.copy_from_slice(d.finalize().as_slice());
    let c_p = Scalar::from_bytes_mod_order_wide(&output);

    let c = c_p * gamma2;

    ((c, c_p, r1, r2, gamma1, gamma2, x, z, a, c), c)
}

pub fn sign_u2(
    stu: (
        Scalar,
        Scalar,
        Scalar,
        Scalar,
        Scalar,
        Scalar,
        RistrettoPoint,
        RistrettoPoint,
        RistrettoPoint,
        Scalar,
    ),
    msg2: (Scalar, Scalar, Scalar),
) -> (Scalar, Scalar, Scalar, Scalar) {
    let (_c, c_p, r1, r2, gamma1, gamma2, _x, _z, _a, _cc) = stu;
    let (s, y, t) = msg2;

    // check y != 0
    // check C == g^t * Z^y
    // check A * X^(cy) == g^s

    let s_p = (gamma1 * gamma2.invert()) * s + r1;
    let y_p = gamma1 * y;
    let t_p = gamma1 * t + r2;

    (c_p, s_p, y_p, t_p)
}

pub fn verify<D>(
    sigma: (Scalar, Scalar, Scalar, Scalar),
    pk: (RistrettoPoint, RistrettoPoint),
    m: &[u8],
) -> bool
where
    D: Digest<OutputSize = U64>,
{
    let (c, s, y, t) = sigma;
    let (x, z) = pk;

    let cc = RistrettoPoint::multiscalar_mul([&t, &y], [&G, &z]);
    let a = RistrettoPoint::multiscalar_mul([&s, &(-c * y)], [&G, &x]);

    let d = D::new()
        .chain_update(a.compress().as_bytes())
        .chain_update(cc.compress().as_bytes())
        .chain_update(m);

    let mut output = [0u8; 64];
    output.copy_from_slice(d.finalize().as_slice());
    let c_p = Scalar::from_bytes_mod_order_wide(&output);

    c == c_p
}

const G: RistrettoPoint = RISTRETTO_BASEPOINT_POINT;

#[cfg(test)]
mod tests {
    use sha3::Sha3_512;

    use super::*;

    #[test]
    fn signing_and_verifying() {
        let (sk, pk) = key_gen();
        let (sts, msg1) = sign_s1(sk, pk);
        let (stu, c) = sign_u1::<Sha3_512>(pk, msg1, b"ok it's fine");
        let msg2 = sign_s2(sts, c);
        let sigma = sign_u2(stu, msg2);
        assert!(verify::<Sha3_512>(sigma, pk, b"ok it's fine"));
    }

    #[test]
    fn bad_message() {
        let (sk, pk) = key_gen();
        let (sts, msg1) = sign_s1(sk, pk);
        let (stu, c) = sign_u1::<Sha3_512>(pk, msg1, b"ok it's fine");
        let msg2 = sign_s2(sts, c);
        let sigma = sign_u2(stu, msg2);
        assert!(!verify::<Sha3_512>(sigma, pk, b"ok it's NOT fine"));
    }
}
