# CloverLeaf_OpenACC

This is the OpenACC version of CloverLeaf version 1.3. 

# Compiling

In most cases one needs to load the appropriate modules (including those for the accelerator in question) and then:

```
make clean; make COMPILER=CRAY
```

or,

```
make clean; make COMPILER=PGI
```
