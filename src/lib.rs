use curve25519_dalek::constants::{RISTRETTO_BASEPOINT_POINT, RISTRETTO_BASEPOINT_TABLE};
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::{IsIdentity, MultiscalarMul};
use digest::consts::U64;
use digest::Digest;
use rand::{CryptoRng, RngCore};

pub struct SecretKey {
    x: Scalar,
    xx: RistrettoPoint,
    zz: RistrettoPoint,
}

impl SecretKey {
    pub fn new<R>(rng: &mut R) -> SecretKey
    where
        R: RngCore + CryptoRng,
    {
        // Algorithm BS_3.KG
        let x = nonzero_scalar(rng);
        let xx = &RISTRETTO_BASEPOINT_TABLE * &x;
        let zz = nonzero_point(rng);

        SecretKey { x, xx, zz }
    }

    pub fn public_key(&self) -> PublicKey {
        PublicKey { xx: self.xx, zz: self.zz }
    }

    pub fn initialize<R>(&self, rng: &mut R) -> (SignerState, SignerInitializeMessage)
    where
        R: RngCore + CryptoRng,
    {
        // Algorithm BS_3.S_1
        let a = nonzero_scalar(rng);
        let t = nonzero_scalar(rng);
        let y = Scalar::random(rng);

        let aa = &RISTRETTO_BASEPOINT_TABLE * &a;
        let cc = RistrettoPoint::multiscalar_mul([&t, &y], [&G, &self.zz]);

        (SignerState { a, y, t, x: self.x }, SignerInitializeMessage { aa, cc })
    }
}

pub struct SignerState {
    a: Scalar,
    y: Scalar,
    t: Scalar,
    x: Scalar,
}

impl SignerState {
    pub fn finalize(self, msg: UserMessage) -> SignerFinalizeMessage {
        // Algorithm BS_3.S_2
        assert_ne!(msg.c, Scalar::zero());

        let s = (msg.c * self.y * self.x) + self.a;

        SignerFinalizeMessage { s, y: self.y, t: self.t }
    }
}

pub struct SignerInitializeMessage {
    aa: RistrettoPoint,
    cc: RistrettoPoint,
}

pub struct SignerFinalizeMessage {
    s: Scalar,
    y: Scalar,
    t: Scalar,
}

pub struct PublicKey {
    xx: RistrettoPoint,
    zz: RistrettoPoint,
}

impl PublicKey {
    pub fn initialize<D, R>(
        &self,
        rng: &mut R,
        msg1: SignerInitializeMessage,
        m: &[u8],
    ) -> (UserState, UserMessage)
    where
        D: Digest<OutputSize = U64>,
        R: RngCore + CryptoRng,
    {
        // Algorithm BS_3.U_1
        let r1 = nonzero_scalar(rng);
        let r2 = nonzero_scalar(rng);
        let gamma1 = Scalar::random(rng);
        let gamma2 = nonzero_scalar(rng);
        // Tessaro & Zhu don't require gamma2 to be non-zero, but inversion of zero is not defined.

        let aa_p =
            RistrettoPoint::multiscalar_mul([&r1, &(gamma1 * gamma2.invert())], [&G, &msg1.aa]);
        let cc_p = RistrettoPoint::multiscalar_mul([&gamma1, &r2], [&msg1.cc, &G]);
        let c_p = hash::<D>(aa_p, cc_p, m);
        let c = c_p * gamma2;

        (
            UserState {
                xx: self.xx,
                zz: self.zz,
                aa: msg1.aa,
                cc: msg1.cc,
                c,
                c_p,
                r1,
                r2,
                gamma1,
                gamma2,
            },
            UserMessage { c },
        )
    }

    pub fn verify<D>(&self, sig: &Signature, m: &[u8]) -> bool
    where
        D: Digest<OutputSize = U64>,
    {
        // Algorithm BS_3.Ver
        if sig.y == Scalar::zero() {
            return false;
        }

        let cc = RistrettoPoint::multiscalar_mul([&sig.t, &sig.y], [&G, &self.zz]);
        let aa = RistrettoPoint::multiscalar_mul([&sig.s, &(-sig.c * sig.y)], [&G, &self.xx]);
        let c_p = hash::<D>(aa, cc, m);

        sig.c == c_p
    }
}

pub struct UserState {
    xx: RistrettoPoint,
    zz: RistrettoPoint,
    aa: RistrettoPoint,
    cc: RistrettoPoint,
    c: Scalar,
    c_p: Scalar,
    r1: Scalar,
    r2: Scalar,
    gamma1: Scalar,
    gamma2: Scalar,
}

impl UserState {
    pub fn finalize(self, msg2: SignerFinalizeMessage) -> Signature {
        // Algorithm BS_3.U_2
        assert_ne!(msg2.y, Scalar::zero());
        assert_eq!(self.cc, RistrettoPoint::multiscalar_mul([&msg2.t, &msg2.y], [&G, &self.zz]));
        assert_eq!(&RISTRETTO_BASEPOINT_TABLE * &msg2.s, self.aa + (self.xx * (self.c * msg2.y)));

        let s_p = (self.gamma1 * self.gamma2.invert()) * msg2.s + self.r1;
        let y_p = self.gamma1 * msg2.y;
        let t_p = self.gamma1 * msg2.t + self.r2;

        Signature { c: self.c_p, s: s_p, y: y_p, t: t_p }
    }
}

pub struct UserMessage {
    c: Scalar,
}

pub struct Signature {
    c: Scalar,
    s: Scalar,
    y: Scalar,
    t: Scalar,
}

#[inline]
fn hash<D>(a: RistrettoPoint, b: RistrettoPoint, m: &[u8]) -> Scalar
where
    D: Digest<OutputSize = U64>,
{
    let d = D::new()
        .chain_update(a.compress().as_bytes())
        .chain_update(b.compress().as_bytes())
        .chain_update(m);

    let mut output = [0u8; 64];
    output.copy_from_slice(d.finalize().as_slice());

    Scalar::from_bytes_mod_order_wide(&output)
}

fn nonzero_scalar<R>(rng: &mut R) -> Scalar
where
    R: RngCore + CryptoRng,
{
    loop {
        let s = Scalar::random(rng);
        if s != Scalar::zero() {
            return s;
        }
    }
}

fn nonzero_point<R>(rng: &mut R) -> RistrettoPoint
where
    R: RngCore + CryptoRng,
{
    loop {
        let p = RistrettoPoint::random(rng);
        if !p.is_identity() {
            return p;
        }
    }
}

const G: RistrettoPoint = RISTRETTO_BASEPOINT_POINT;

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_chacha::ChaChaRng;
    use sha3::Sha3_512;

    use super::*;

    #[test]
    fn signing_and_verifying() {
        let mut rng = ChaChaRng::seed_from_u64(100);
        let msg = b"this is a secret";
        let sk = SecretKey::new(&mut rng);
        let (s1, sm1) = sk.initialize(&mut rng);
        let pk = sk.public_key();
        let (u1, um1) = pk.initialize::<Sha3_512, _>(&mut rng, sm1, msg);
        let sm2 = s1.finalize(um1);
        let sig = u1.finalize(sm2);

        assert!(pk.verify::<Sha3_512>(&sig, msg));
        assert!(!pk.verify::<Sha3_512>(&sig, b"some other guy"));
    }
}
