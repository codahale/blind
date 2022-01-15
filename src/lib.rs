use std::marker::PhantomData;

use curve25519_dalek::constants::{RISTRETTO_BASEPOINT_POINT, RISTRETTO_BASEPOINT_TABLE};
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::traits::MultiscalarMul;
use digest::consts::U64;
use digest::Digest;
use rand::{CryptoRng, RngCore};

pub struct SecretKey<D: Digest<OutputSize = U64>> {
    x: Scalar,
    xx: RistrettoPoint,
    zz: RistrettoPoint,
    _hash: PhantomData<D>,
}

impl<D: Digest<OutputSize = U64>> SecretKey<D> {
    pub fn new<R: RngCore + CryptoRng>(rng: &mut R) -> SecretKey<D> {
        let x = Scalar::random(rng);
        let xx = &RISTRETTO_BASEPOINT_TABLE * &x;
        let zz = RistrettoPoint::random(rng);

        SecretKey { x, xx, zz, _hash: Default::default() }
    }

    pub fn public_key(&self) -> PublicKey<D> {
        PublicKey { xx: self.xx, zz: self.zz, _hash: Default::default() }
    }

    pub fn initialize<R: RngCore + CryptoRng>(
        &self,
        rng: &mut R,
    ) -> (SignerState<D>, SignerMessage1<D>) {
        let a = Scalar::random(rng);
        let t = Scalar::random(rng);
        let y = Scalar::random(rng);

        let aa = &RISTRETTO_BASEPOINT_TABLE * &a;
        let cc = RistrettoPoint::multiscalar_mul([&t, &y], [&G, &self.zz]);

        (
            SignerState { a, y, t, x: self.x, _hash: Default::default() },
            SignerMessage1 { aa, cc, _hash: Default::default() },
        )
    }
}

pub struct SignerState<D: Digest<OutputSize = U64>> {
    a: Scalar,
    y: Scalar,
    t: Scalar,
    x: Scalar,
    _hash: PhantomData<D>,
}

impl<D: Digest<OutputSize = U64>> SignerState<D> {
    pub fn finalize(self, msg: UserMessage<D>) -> SignerMessage2<D> {
        let s = (msg.c * self.y * self.x) + self.a;
        SignerMessage2 { s, y: self.y, t: self.t, _hash: Default::default() }
    }
}

pub struct SignerMessage1<D: Digest<OutputSize = U64>> {
    aa: RistrettoPoint,
    cc: RistrettoPoint,
    _hash: PhantomData<D>,
}

pub struct SignerMessage2<D: Digest<OutputSize = U64>> {
    s: Scalar,
    y: Scalar,
    t: Scalar,
    _hash: PhantomData<D>,
}

pub struct PublicKey<D: Digest<OutputSize = U64>> {
    xx: RistrettoPoint,
    zz: RistrettoPoint,
    _hash: PhantomData<D>,
}

impl<D: Digest<OutputSize = U64>> PublicKey<D> {
    pub fn initialize<R: RngCore + CryptoRng>(
        &self,
        rng: &mut R,
        msg1: SignerMessage1<D>,
        m: &[u8],
    ) -> (UserState<D>, UserMessage<D>) {
        let r1 = Scalar::random(rng);
        let r2 = Scalar::random(rng);
        let gamma1 = Scalar::random(rng);
        let gamma2 = Scalar::random(rng);

        let aa =
            RistrettoPoint::multiscalar_mul([&r1, &(gamma1 * gamma2.invert())], [&G, &msg1.aa]);
        let cc = RistrettoPoint::multiscalar_mul([&gamma1, &r2], [&msg1.cc, &G]);

        let d = D::new()
            .chain_update(aa.compress().as_bytes())
            .chain_update(cc.compress().as_bytes())
            .chain_update(m);

        let mut output = [0u8; 64];
        output.copy_from_slice(d.finalize().as_slice());
        let c_p = Scalar::from_bytes_mod_order_wide(&output);

        let c = c_p * gamma2;

        (
            UserState { c_p, r1, r2, gamma1, gamma2, _hash: Default::default() },
            UserMessage { c, _hash: Default::default() },
        )
    }

    pub fn verify(&self, sig: &Signature<D>, m: &[u8]) -> bool {
        let cc = RistrettoPoint::multiscalar_mul([&sig.t, &sig.y], [&G, &self.zz]);
        let a = RistrettoPoint::multiscalar_mul([&sig.s, &(-sig.c * sig.y)], [&G, &self.xx]);

        let d = D::new()
            .chain_update(a.compress().as_bytes())
            .chain_update(cc.compress().as_bytes())
            .chain_update(m);

        let mut output = [0u8; 64];
        output.copy_from_slice(d.finalize().as_slice());
        let c_p = Scalar::from_bytes_mod_order_wide(&output);

        sig.c == c_p
    }
}

pub struct UserState<D: Digest<OutputSize = U64>> {
    c_p: Scalar,
    r1: Scalar,
    r2: Scalar,
    gamma1: Scalar,
    gamma2: Scalar,
    _hash: PhantomData<D>,
}

impl<D: Digest<OutputSize = U64>> UserState<D> {
    pub fn finalize(self, msg2: SignerMessage2<D>) -> Signature<D> {
        let s_p = (self.gamma1 * self.gamma2.invert()) * msg2.s + self.r1;
        let y_p = self.gamma1 * msg2.y;
        let t_p = self.gamma1 * msg2.t + self.r2;

        Signature { c: self.c_p, s: s_p, y: y_p, t: t_p, _hash: Default::default() }
    }
}

pub struct UserMessage<D: Digest<OutputSize = U64>> {
    c: Scalar,
    _hash: PhantomData<D>,
}

pub struct Signature<D: Digest<OutputSize = U64>> {
    c: Scalar,
    s: Scalar,
    y: Scalar,
    t: Scalar,
    _hash: PhantomData<D>,
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

        let sk = SecretKey::<Sha3_512>::new(&mut rng);

        let (s1, sm1) = sk.initialize(&mut rng);

        let pk = sk.public_key();

        let (u1, um1) = pk.initialize(&mut rng, sm1, msg);

        let sm2 = s1.finalize(um1);

        let sig = u1.finalize(sm2);

        assert!(pk.verify(&sig, msg));
        assert!(!pk.verify(&sig, b"some other guy"));
    }
}
