#   mmKyber-artifact

A multi-message, multi-recipient PKE/KEM enables the batch encryption of
individual messages for multiple independent recipients in a single operation,
reducing costs, particularly bandwidth, compared to the trivial solution of
encrypting each message individually.

What's where:
```
mmKyber-artifact
├── mmKyber-c         # Plain C Implementation of mmKyber. Used for benchmarks.
├── mmKyber-py        # Python model for mmKyber, compatible with the C code.
├── mmKyber-pkzk      # ZK proofs of mmKyber public key validity with LaZer.
├── refKyber-c        # Kyber reference and "apples-to-apples" benchmarks.
└── README.md         # this file
```

