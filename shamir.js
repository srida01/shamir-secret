// shamir.js — Shamir's Secret Sharing Scheme (GF prime field, BigInt)

export const PRIME = 208351617316091241234326746312124448251235562226470491514186331217050270460481n;

export function mod(a, m) {
  return ((a % m) + m) % m;
}

export function modPow(base, exp, m) {
  let result = 1n;
  base = mod(base, m);
  while (exp > 0n) {
    if (exp % 2n === 1n) result = result * base % m;
    exp = exp / 2n;
    base = base * base % m;
  }
  return result;
}

export function modInv(a, m) {
  // Fermat's little theorem: a^(p-2) mod p
  return modPow(a, m - 2n, m);
}

export function secretToInt(str) {
  let n = 0n;
  for (let i = 0; i < str.length; i++) {
    n = n * 256n + BigInt(str.charCodeAt(i));
  }
  return n;
}

export function intToSecret(n) {
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

// Shorten a big BigInt for display
function sh(b, head = 12, tail = 6) {
  const s = b.toString();
  return s;
}

/**
 * Split a secret string into n shares requiring k to reconstruct.
 * @returns {{ shares: {x:bigint,y:bigint}[], steps: object[] }}
 */
export function split(secret, n, k) {
  if (!secret) throw new Error('--secret is required');
  if (k < 2)   throw new Error('threshold k must be ≥ 2');
  if (k > n)   throw new Error('threshold k cannot exceed shares n');
  if (n > 20)  throw new Error('max 20 shares supported');

  const steps = [];

  // ── Step 1: Encode ────────────────────────────────────────────────────────
  const s = secretToInt(secret);
  const charBytes = [...secret].map(c => `'${c}'=0x${c.charCodeAt(0).toString(16)}`).join(', ');
  steps.push({
    title: 'Step 1 — Encode secret as integer',
    color: 'accent2',
    lines: [
      { type: 'formula', text: 'secret  →  UTF-8 bytes  →  base-256 big integer' },
      { type: 'sub',     text: 'Each char contributes 8 bits: n = n×256 + charCode(c)' },
      { type: 'data',    label: 'chars',  text: charBytes },
      { type: 'data',    label: 's (int)', text: sh(s) },
      { type: 'note',    text: 's becomes the constant term a₀ of our secret polynomial f(x).' },
    ]
  });

  // ── Step 2: Build polynomial ───────────────────────────────────────────────
  const coeffs = [s];
  for (let i = 1; i < k; i++) coeffs.push(1n + randomBigIntBelow(PRIME - 2n));

  const polyStr = coeffs.map((c, i) => {
    if (i === 0) return 'a₀';
    if (i === 1) return 'a₁x';
    return `a${i}xᵃ`.replace('ᵃ', i > 9 ? i : ['⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹'][i]);
  }).join(' + ');

  steps.push({
    title: `Step 2 — Construct random degree-${k - 1} polynomial`,
    color: 'accent',
    lines: [
      { type: 'formula', text: `f(x) = ${polyStr}  (mod p)` },
      { type: 'sub',     text: `Degree k-1 = ${k - 1}. Any k points uniquely determine it; fewer points reveal nothing about f(0).` },
      { type: 'data',    label: 'a₀ = s',    text: sh(s) + '  ← secret (constant term)' },
      ...coeffs.slice(1).map((c, i) => ({
        type: 'data', label: `a${i + 1}`, text: sh(c) + '  ← random in [1, p-1]'
      })),
      { type: 'data',    label: 'prime p',   text: sh(PRIME, 14, 8) },
    ]
  });

  // ── Step 3: Evaluate ───────────────────────────────────────────────────────
  const shares = [];
  const evalLines = [];
  for (let x = 1; x <= n; x++) {
    let y = 0n;
    const xb = BigInt(x);
    const termParts = coeffs.map((c, i) => {
      const val = mod(c * modPow(xb, BigInt(i), PRIME), PRIME);
      y = mod(y + val, PRIME);
      return `a${i}·${x}^${i}`;
    });
    shares.push({ x: xb, y });
    evalLines.push({
      type: 'data',
      label: `f(${x})`,
      text: `${termParts.join(' + ')} mod p  =  ${sh(y)}`
    });
  }

  steps.push({
    title: `Step 3 — Evaluate f(x) at x = 1 … ${n}  →  shares`,
    color: 'ok',
    lines: [
      { type: 'formula', text: 'share_i = ( i,  f(i) mod p )  for i = 1 … n' },
      { type: 'sub',     text: 'Each party receives one (x, y) pair. The x value is public; y is the share secret.' },
      ...evalLines,
      { type: 'note',    text: `Distribute these ${n} shares. Any ${k} can reconstruct f(0) = s.` },
    ]
  });

  return { shares, steps };
}

/**
 * Reconstruct the secret from an array of shares via Lagrange interpolation.
 * @returns {{ secret: string, steps: object[] }}
 */
export function reconstruct(shares) {
  if (!shares || shares.length === 0) throw new Error('no shares provided');

  const steps = [];
  const k = shares.length;

  // ── Step 1: Setup ─────────────────────────────────────────────────────────
  steps.push({
    title: 'Step 1 — Setup: Lagrange interpolation at x = 0',
    color: 'accent2',
    lines: [
      { type: 'formula', text: 'f(0) = Σᵢ  yᵢ · Lᵢ(0)  (mod p)' },
      { type: 'sub',     text: `Using ${k} shares. Each basis polynomial Lᵢ(0) is 1 at xᵢ and 0 at all other xⱼ.` },
      { type: 'formula', text: 'Lᵢ(0) = Π_{j≠i} (0 − xⱼ) / (xᵢ − xⱼ)  mod p' },
      { type: 'sub',     text: 'Division mod p uses modular inverse via Fermat\'s little theorem: a⁻¹ ≡ a^(p−2) mod p' },
      ...shares.map((s, i) => ({
        type: 'data', label: `share[${i + 1}]`, text: `x = ${s.x},  y = ${sh(s.y)}`
      })),
    ]
  });

  // ── Step 2: Basis polynomials ─────────────────────────────────────────────
  const bases = [];
  const basisLines = [];
  for (let i = 0; i < k; i++) {
    let num = 1n, den = 1n;
    const numParts = [], denParts = [];
    for (let j = 0; j < k; j++) {
      if (i === j) continue;
      numParts.push(`(−${shares[j].x})`);
      denParts.push(`(${shares[i].x}−${shares[j].x})`);
      num = mod(num * mod(-shares[j].x, PRIME), PRIME);
      den = mod(den * mod(shares[i].x - shares[j].x, PRIME), PRIME);
    }
    const L = mod(num * modInv(den, PRIME), PRIME);
    bases.push(L);
    basisLines.push({
      type: 'data',
      label: `L${i}(0)`,
      text: `num=${numParts.join('·')}  /  den=${denParts.join('·')}  =  ${sh(L)}`
    });
  }

  steps.push({
    title: 'Step 2 — Compute each Lagrange basis Lᵢ(0)',
    color: 'accent',
    lines: [
      { type: 'formula', text: 'Lᵢ(0) = Π_{j≠i} (−xⱼ) · [ Π_{j≠i} (xᵢ−xⱼ) ]⁻¹  mod p' },
      ...basisLines,
    ]
  });

  // ── Step 3: Weighted sum ──────────────────────────────────────────────────
  const termLines = [];
  let secret = 0n;
  for (let i = 0; i < k; i++) {
    const term = mod(shares[i].y * bases[i], PRIME);
    termLines.push({
      type: 'data',
      label: `y${i}·L${i}(0)`,
      text: `${sh(shares[i].y)} · ${sh(bases[i])} mod p  =  ${sh(term)}`
    });
    secret = mod(secret + term, PRIME);
  }

  steps.push({
    title: 'Step 3 — Weighted sum → f(0) = secret integer',
    color: 'ok',
    lines: [
      { type: 'formula', text: 'f(0) = ( Σᵢ yᵢ·Lᵢ(0) ) mod p' },
      ...termLines,
      { type: 'data',    label: 'f(0)',   text: sh(secret) + '  ← encoded secret' },
    ]
  });

  // ── Step 4: Decode ────────────────────────────────────────────────────────
  const recovered = intToSecret(secret);
  steps.push({
    title: 'Step 4 — Decode integer back to string',
    color: 'accent2',
    lines: [
      { type: 'formula', text: 'integer  →  base-256 bytes  →  UTF-8 string' },
      { type: 'sub',     text: 'Repeat: charCode = int mod 256, int = int ÷ 256, until 0' },
      { type: 'data',    label: 'integer', text: sh(secret) },
      { type: 'data',    label: 'secret',  text: `"${recovered}"` },
    ]
  });

  return { secret: recovered, steps };
}