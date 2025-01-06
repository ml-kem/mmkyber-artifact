#   mmKyber-artifact

A multi-message multi-recipient PKE/KEM enables the batch encryption of
multiple messages (as a message vector) for multiple independent recipients
in a single operation, significantly reducing costs, particularly bandwidth,
compared to the trivial solution of encrypting each message individually.

What's here:
```
mmKyber-artifact
├── mmKyber-c       # Plain C Implementation of mmKyber. Used for benchmarks.
├── mmKyber-py      # Python model for mmKyber, compatible with the C code.
├── mmKyber-pkzk    # ZK proofs of public key validity with LaZer.
├── pr-fail-dec     # Computation of decryption/decaps failure probabilities.
├── refKyber-c      # Plain C Kyber ref code for "apples-to-apples" benchmarks.
└── README.md       # this file
```

