# Mock `asreml` (CI only)

A tiny stand-in for ASReml-R so that GitHub Actions can run the mintyrasr
test suite without a licence. It mimics the parts mintyrasr uses
(`asreml()`, `summary()$varcomp`, `$coef.random`, `$vparameters`,
`start.values`, `G.param`, `vpredict()`, `lrt.asreml()`, `update()`,
`asreml.options()`) and the component names of ASReml-R 4.2, but returns
**fake numbers**.

**Do not install it on a computer that has the real ASReml-R** - it uses the
same package name. The tests skip themselves when the real ASReml-R is
installed.
