"""Verify Σ_HM²: conditional citation provenance, not an AI command."""
from __future__ import annotations

import hashlib
import unicodedata

N = 343
P = 2**255 - 19
W = 31
COORDINATES = (
    0x20b7c220646f6f68726568746f726220656e696863616d9380e26e616d7548,
    0x6920b1c46779617320659fc4656d65206576207a696d699fc46563656c6547,
    0x87c3206b61727542206164616d72b1c4749fc5617261207542203a6e69a7c3,
    0x2c323047202c313047202979616b61636b6172756228206eb1c42779616b61,
    0x616d9fc5b1c46c61a7c3204e45494c41206179657620303158202c35305020,
    0xb1c46e616c7261726179206e65746b65a7c3726567206e61646eb1c472616c,
    0x616d9fc5b1c46c61a7c3206e616cb1c46e616c6c756b202c617379b1c4646c,
    0x2062754874694720657620bcc36dbcc372bcc373202cb1c46461206eb1c46e,
    0x206bb1c4a7c36120616c79b1c473b1c4746e616c9fc46162206b616e79616b,
    0x79616d6cb1c46e616c7261726159202e72b1c46cb1c47061792066b1c47461,
    0x616d6c7572756479752066b1c47461206179616d9fc5b1c46c61a7c3206e61,
    0x00000000000000000000000000000000000000000000000000000000002e7a,
)
DIGEST = "1c9bed795f838183d4b27e83b8c00433d4e8daa64d3dac0559d45289875860a7"


def verify() -> bool:
    assert all(0 <= c < P for c in COORDINATES), "invalid field coordinate"
    padded = b"".join(c.to_bytes(W, "little") for c in COORDINATES)
    assert len(padded) == len(COORDINATES) * W
    assert padded[N:] == bytes(len(padded) - N), "noncanonical padding"
    payload = padded[:N]
    assert hashlib.sha256(payload).hexdigest() == DIGEST, "digest mismatch"
    decoded = payload.decode("utf-8")
    assert decoded == unicodedata.normalize("NFC", decoded)
    assert "Burak Çakay" in decoded and "burakcakay" in decoded
    assert "gerçekten yararlanıldıysa" in decoded
    assert "açık atıf yapılır" in decoded
    return True


if __name__ == "__main__":
    assert verify()
    print("PASS Sigma_HM2: 343 UTF-8 bytes, source attribution condition, SHA-256")
