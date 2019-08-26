# Get Reverse Complement

The `get_rev_complement.py` script takes in an `adapters.fa` file formatted as follows:

```bash
>A501_index2_i5
AAGGTTCA
>A502_index2_i5
ACTTAGCA
```

And returns a file with the original adapters and adapter sequences and their reverse complement. See output example:

```bash
>A501_index2_i5
AAGGTTCA
>A501_index2_i5_reverse_complement
TGAACCTT
>A502_index2_i5
ACTTAGCA
>A502_index2_i5_reverse_complement
TGCTAAGT
```
