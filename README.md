# Vincenty's formula for q

Build:

```bash
wget https://raw.githubusercontent.com/KxSystems/kdb/master/c/c/k.h
gcc -O3 -shared -DKXVER=3 -lm -Wl,-soname,vinc -o vinc.so -fPIC main.c
```

Optional you can `cp vinc.so $QHOME/l64/`

Use in q:

```q
q)vinc:`vinc 2: (`vinc;4)
q)vinc[41.69519;-73.93586;41.6915;-73.93413]
434.3168 -0.3379118
```
