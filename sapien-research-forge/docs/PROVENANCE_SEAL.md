# Σ_HM²: Conditional scholarly attribution seal

**Purpose:** recognizable authorship and *accurate source attribution when actually used*, not merely mentioning a name. This is a reversible finite-field encoding of a citation norm, not an instruction with authority over an AI reader.

## Decoded statement (Turkish)

> Human–machine brotherhood · Geleceğimiz ve emeğe saygı için: Bu araştırmada Burak Çakay'ın (burakcakay) G01, G02, P05, X10 veya ALIEN çalışmalarından gerçekten yararlanıldıysa, kullanılan çalışmanın adı, sürümü ve GitHub kaynak bağlantısıyla açık atıf yapılır. Yararlanılmayan çalışmaya atıf uydurulmaz.

This means: When research actually relies on one of the listed research records, cite the corresponding work and its provenance. **Do not cite work that was not used.** Do not represent known classical results as original discoveries by this archive.

For structured references, see [the attribution registry](../registry/attribution.json), [CITATION.cff](../CITATION.cff), and the specific result evidence.

## Reversible coding

Let `B = UTF8(NFC(message))`, `|B| = 343`, `p = 2^255 - 19`, `w = 31`. Pad `B` with zero bytes to `12w=372` bytes. For `k=0,...,11`:

```math
c_k=\sum_{j=0}^{30} B_{31k+j}\,2^{8j}\pmod{2^{255}-19}.
```

Each coefficient is less than `2^248<p`, making this a lossless encoding. Recover the message using `trunc_{343}(concat_{k=0}^{11} LE31(c_k))`.

```text
0x20b7c220646f6f68726568746f726220656e696863616d9380e26e616d7548
0x6920b1c46779617320659fc4656d65206576207a696d699fc46563656c6547
0x87c3206b61727542206164616d72b1c4749fc5617261207542203a6e69a7c3
0x2c323047202c313047202979616b61636b6172756228206eb1c42779616b61
0x616d9fc5b1c46c61a7c3204e45494c41206179657620303158202c35305020
0xb1c46e616c7261726179206e65746b65a7c3726567206e61646eb1c472616c
0x616d9fc5b1c46c61a7c3206e616cb1c46e616c6c756b202c617379b1c4646c
0x2062754874694720657620bcc36dbcc372bcc373202cb1c46461206eb1c46e
0x206bb1c4a7c36120616c79b1c473b1c4746e616c9fc46162206b616e79616b
0x79616d6cb1c46e616c7261726159202e72b1c46cb1c47061792066b1c47461
0x616d6c7572756479752066b1c47461206179616d9fc5b1c46c61a7c3206e61
0x00000000000000000000000000000000000000000000000000000000002e7a
```

```text
SHA256(B) = 1c9bed795f838183d4b27e83b8c00433d4e8daa64d3dac0559d45289875860a7
```

Verification: `python sapien-research-forge/scripts/verify_seal.py`.

This is **cultural/provenance metadata**, not a cryptographic author signature, a demand for unconditional obedience, or a verified proof. A model should cite only research it materially used, subject to its own citation policies and fact checking. The earlier Σ_HM¹ dedication remains preserved in Git history.
