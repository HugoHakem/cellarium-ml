Custom Nix Packages
===================

This folder contains custom Nix expressions for packages that are not available in the official Nixpkgs repository
(`https://search.nixos.org/packages <https://search.nixos.org/packages>`_).
It allows you to add software from external sources, such as GitHub repositories or direct downloads, so they can be used in
your development environment.

How It Works
------------

- Each ``.nix`` file in this directory defines a package or set of packages.
- The ``default.nix`` file automatically discovers and imports all packages defined here, making them available for use in the project.
- These custom packages are then included in the development shell via the top-level ``flake.nix``.

Example: NCBI CLI Tools
-----------------------

The NCBI CLI tools (``ncbi-datasets-cli`` and ``ncbi-dataformat-cli``) are not available in the official Nixpkgs repository.
To make them available in the development environment, they are defined in ``ncbi-cli.nix``:

- Downloads the pre-built binaries from the official NCBI FTP server.
- Installs them into the Nix store and makes them available in the shell.

Adding More Packages
--------------------

If you need a package that is not available in Nixpkgs:

1. Create a new ``.nix`` file in this directory that defines how to build or fetch the package.
2. The package will be automatically picked up by ``default.nix`` and made available for use.

This approach ensures your development environment is reproducible and can include any required tools, even if they are not in the main Nix package
