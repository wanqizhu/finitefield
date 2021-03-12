# Finite Field for Error Correcting Codes

This project implements a generic FiniteField class, supporting Galois fields
of any (prime powered) size. As an example for how one might use this class for
Error Correcting Codes (ECC), we provide a sample implementation of the Reed
Solomon code, using the Berlekamp Welch algorithm as an encoder.

While correct functionality is necessary, the goal of this project is to build
general purpose, well-decomposed classes that are easy to use and modify. In
particular, the Finite Field class has no knowledge or dependency on how it is
to be used. One can easily substitute in other ECC implementations or other
mathematical uses on top of the Finite Field module.

Lastly, this project has zero dependencies. While using packages such as
`numpy` would simplify some implementation, we designed the classes so that
they are completely self-contained. This allows greater customization for users
who want to completely understand or modify the internal implementations.


# Design Overview

We provide the following classes:

### [FiniteField](FiniteField.py)
Construct a finite field of a given order. For prime ordered fields, the field
elements are simply the integers `0` to `p-1`, with addition and multiplication
mod p. For fields of order `p^m` for `m > 1`, construction requires passing in
a primitive polynomial `g`, and the field elements are the polynomial
remainders `F_p[x] mod g`.

```
F5 = FiniteField(5)                 # order 5 
F49 = FiniteField(7, 2, [1, 0, 1])  # order 7^2 = 49, with primitive polynomial
                                    # 1 + x^2 (1 + 0 x + 1 x^2)
```

After construction, one can use an instance of this class like a function to
create field elements (see the FiniteFieldElem class below):

```
print(F5(3))  # just the number 3, in the context of the integers mod 5
assert F5(2) + F5(4) == F5(1)  # 2 + 4 == 1 mod 5

# each element of F_{49} is a remainder polynomial of degree < 2 in F_7[x]
a = F49([1, 3])   # the element 1 + 3x
b = F49([3, 4])
print(a * b)      # should get [F49] (5, 6), since (1 + 3x) * (3 + 4x) mod (1 + x^2)
                  # == 3 + 13x + 12x^2
                  # == -9 + 13x
                  # == 5 + 6x  (remember each coefficient is taken mod 7)
```

For a mathematical walkthrough of finite fields of nonprime order, see [these
lecture notes](http://web.stanford.edu/~marykw/classes/CS250_W19/readings/Forney_Introduction_to_Finite_Fields.pdf).

By default, on field creation a linear search is performed to find a primitive
element of this field. The field object returned then keeps track of an
internal log table, which allows `O(1)` multiplication and division. This
behavior can be turned off by passing in `find_primitive_elem = False` to the
constructor.

### [FiniteFieldElem](FiniteField.py#L200)
Internal class representing a field element. Rather than creating an instance
of this class directly, it is recommended to use the field as a function call,
as shown in the example above.

This class overwrites the arithmetic operators `+, -, *, /`. These operators
are defined over two elements of the same field. You can also directly operate with an integer, in which case the integer is interpretted as a field element with all higher ordered terms set to `0`.


### [ReedSolomon](ReedSolomon.py)

This class gives a working example for how one might use the FiniteField class.
In particular, we implement the Reed Solomon code, an Error Correcting Code
that achives the Singleton Bound and is widely used in practice.

Note that this implementation requires no knowledge about the internal workings
of the FiniteField or FiniteFieldElem class. The only requirement is that the
message to be encoded or decoded are given as lists of field elements, and that
those field elements support arithmetic operators. This demonstrates the
decomposability and general-purposeness of our design.

-----

In addition, we provide a library of general-purpose helper math functions,
such as Gaussian Elimination and Polynomial Interpolation, in
[utils.py](utils.py). These functions operate over any data type that supports
basic arithmetic `(+, -, *, /)`, so one can use them over the reals, over
FiniteFieldElems, and over any alternative field implementations.

# Testing

Unittests are defined in [test.py](test.py). Run all tests with

`python -m unittest`
