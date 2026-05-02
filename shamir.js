// shamir.js — Shamir's Secret Sharing Scheme (GF prime field, BigInt)

const PRIME = 208351617316091241234326746312124448251235562226470491514186331217050270460481n;

function mod(a, m) {
  return ((a % m) + m) % m;
}

function modPow(base, exp, m) {
  let result = 1n;
  base = mod(base, m);
  while (exp > 0n) {
    if (exp % 2n === 1n) result = result * base % m;
    exp = exp / 2n;
    base = base * base % m;
  }
  return result;
}

function modInv(a, m) {
  // Fermat's little theorem: a^(p-2) mod p
  return modPow(a, m - 2n, m);
}

function secretToInt(str) {
  let n = 0n;
  for (let i = 0; i < str.length; i++) {
    n = n * 256n + BigInt(str.charCodeAt(i));
  }
  return n;
}

function intToSecret(n) {
  let str = '';
  while (n > 0n) {
    str = String.fromCharCode(Number(n % 256n)) + str;
    n = n / 256n;
  }
  return str;
}

function randomBigIntBelow(max) {
  const byteLen = Math.ceil(max.toString(16).length / 2) + 4;
  const bytes = new Uint8Array(byteLen);
  crypto.getRandomValues(bytes);
  let n = 0n;
  for (const b of bytes) n = n * 256n + BigInt(b);
  return n % max;
}

/**
 * Split a secret string into n shares requiring k to reconstruct.
 * @param {string} secret
 * @param {number} n - total shares
 * @param {number} k - threshold
 * @returns {{ x: bigint, y: bigint }[]}
 */
export function split(secret, n, k) {
  if (!secret) throw new Error('--secret is required');
  if (k < 2) throw new Error('threshold k must be ≥ 2');
  if (k > n) throw new Error('threshold k cannot exceed shares n');
  if (n > 20) throw new Error('max 20 shares supported');

  const s = secretToInt(secret);
  // Random polynomial: f(x) = s + a1*x + ... + a(k-1)*x^(k-1)
  const coeffs = [s];
  for (let i = 1; i < k; i++) {
    coeffs.push(1n + randomBigIntBelow(PRIME - 2n));
  }

  const shares = [];
  for (let x = 1; x <= n; x++) {
    let y = 0n;
    const xb = BigInt(x);
    for (let i = 0; i < k; i++) {
      y = mod(y + coeffs[i] * modPow(xb, BigInt(i), PRIME), PRIME);
    }
    shares.push({ x: xb, y });
  }
  return shares;
}

/**
 * Reconstruct the secret from an array of shares via Lagrange interpolation.
 * @param {{ x: bigint, y: bigint }[]} shares
 * @returns {string}
 */
export function reconstruct(shares) {
  if (!shares || shares.length === 0) throw new Error('no shares provided');

  let secret = 0n;
  for (let i = 0; i < shares.length; i++) {
    let num = shares[i].y;
    let den = 1n;
    for (let j = 0; j < shares.length; j++) {
      if (i === j) continue;
      num = mod(num * mod(-shares[j].x, PRIME), PRIME);
      den = mod(den * mod(shares[i].x - shares[j].x, PRIME), PRIME);
    }
    secret = mod(secret + num * modInv(den, PRIME), PRIME);
  }
  return intToSecret(secret);
}