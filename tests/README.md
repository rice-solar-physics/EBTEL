# Testing with Python

This directory contains a number of tests that can be run with the
[`pytest`](https://docs.pytest.org/en/stable/) testing framework in Python.
To install the needed dependencies,

```shell
$ pip install pytest numpy hissw matplotlib
```

You will also need to properly configure [`hissw`](https://wtbarnes.github.io/hissw/) for your particular IDL installation. Then, to run the tests,

```shell
$ pytest --ebtel_idl_path=/path/to/EBTEL
```

where `/path/to/EBTEL` is the path to the root of the EBTEL repository on your machine.
Alternatively, you can also plot the results from the tests to visually confirm
that everything is working as expected,

```shell
$ pytest --ebtel_idl_path=/path/to/EBTEL --show_plots
```

Note that these tests do not necessarily test the correctness of the code, but simply
test that the code runs without crashing for several different input configurations.
