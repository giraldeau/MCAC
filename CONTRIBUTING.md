# Contributing

First off, thank you for considering contributing to MCAC.

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

Following these guidelines helps to communicate that you respect
the time of the developers managing and developing this project.
In return, they should reciprocate that respect in addressing your issue, assessing changes,
and helping you finalize your pull requests.

For any questions not concerning directly MCAC, 
read the docs, ask Google, look on Stack Overflow (in this order).

If your problem is indeed MCAC specific,
ask one of the maintainers by writing an [issue on GitLab](https://gitlab.coria-cfd.fr/MCAC/MCAC/issues).

The current stable version is in the branch `master`  
The branch `next` contains beta developments that will regularly be merged with `master` for each new version.   
All other future developments are on their respective branches and will be merged into `next`

## Types of Contributions

You can contribute in many ways:

### Report Bugs

Report bugs directly our [issues list on GitLab](https://gitlab.coria-cfd.fr/MCAC/MCAC/issues).

Please include,
 * the version of MCAC
 * detailed steps to reproduce the bug.

If you don't have steps to reproduce the bug,
just note your observations in as much detail as you can.

### Suggest features or enhancement

If you find yourself wishing for a feature that doesn't exist in MCAC,
you are probably not alone. There are bound to be others out there with similar needs.

Any code or part of a script that you expect to be used multiple times
is relevelant for being included in this library.  

Open an issue on our [issues list on GitLab](https://gitlab.coria-cfd.fr/MCAC/MCAC/issues): 

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.

### Fix bugs

Look through the [GitLab issues](https://gitlab.coria-cfd.fr/MCAC/MCAC/issues) for bugs.

### Implement features or enhancement

Look through the [GitLab issues](https://gitlab.coria-cfd.fr/MCAC/MCAC/issues) for missing features.
Please do not combine multiple feature enhancements into a single pull request.

## Setting Up the Code for Local Development

1. Clone the project locally.  
   ```shell
   git clone git@gitlab.coria-cfd.fr:MCAC/MCAC.git
   ```

2. Compile your local copy as usual.  
   ```shell 
   cd MCAC/src  
   qmake MCAC.pro
   make
   cd ..
   pip install -e .
   ```

3. Create a branch for local development.  
   ```shell
   git checkout -b name-of-your-bugfix-or-feature
   ```
    
4. Make a lot's of little commits while you work
   ```shell
   git add ...  
   git commit -m "Your detailed description of your changes"  
   git push origin name-of-your-bugfix-or-feature
   ```  

5. Document your work using docstrings

6. When you are fully satisfied
 - merge the current state of the `next` branch  
   ```shell
   git merge next
   ```  
 - make a [push request](https://gitlab.coria-cfd.fr/MCAC/MCAC/branches)  
   which will alert the main contributers and they will finish to merge your developments.


[comment]: <> (5. Checks that the tests are still working)
[comment]: <> (7. If this is a new feature, write one or multiple tests cases)

## Contributor Guidelines

### Coding Standards

* English language 

[comment]: <> (### Testing)

## Core Committer Guide
