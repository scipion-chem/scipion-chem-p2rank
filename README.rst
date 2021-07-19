=======================
P2Rank plugin
=======================

This is a **Scipion** plugin that offers different **p2rank tools**.
These tools will make it possible to carry out different functions for predicting protein pockets

Therefore, this plugin allows to use programs from the p2rank software suite
within the Scipion framework.

==========================
Install this plugin
==========================

You will need to use `Scipion3 <https://scipion-em.github.io/docs/docs/scipion
-modes/how-to-install.html>`_ to run these protocols.


1. **Install the plugin in Scipion**

P2Rank is installed automatically by scipion.

- **Install the stable version (Not available yet)**

    Through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

    or

.. code-block::

    scipion3 installp -p scipion-chem-p2rank


- **Developer's version**

    1. Download repository:

    .. code-block::

        git clone https://github.com/scipion-chem/scipion-chem-p2rank.git

    2. Install:

    .. code-block::

        scipion3 installp -p path_to_scipion-chem-p2rank --devel


