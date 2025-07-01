Custom Packages Directory
========================

This folder is intended for defining custom Nix packages that are not available in the `Nix Packages collection <https://search.nixos.org/packages>`_.

When to Use
-----------

- **Unavailable Packages:** If you need a package or tool (for example, from a GitHub repository) that cannot be found in the official Nixpkgs repository, you can add it here.
- **Custom Builds:** If you want to override the build process or patch a package in a way that is not possible through overlays or upstream Nixpkgs.

How It Works
------------

- Each ``.nix`` file in this directory (except ``default.nix``) should define one or more packages as Nix attributes.
- The ``default.nix`` file automatically imports all other ``.nix`` files in this directory and aggregates their outputs.
- The resulting set of packages is imported in the main ``flake.nix`` and can be added to your development shell or used elsewhere in your project.

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

Further Example
---------------

To add a custom package, create a new file (e.g., ``mytool.nix``) in this directory:

.. code-block:: nix

   self: super:

   {
     mytool = super.stdenv.mkDerivation {
       pname = "mytool";
       version = "1.0.0";
       src = super.fetchFromGitHub {
         owner = "username";
         repo = "mytool";
         rev = "commit-or-tag";
         sha256 = "sha256-...";
       };
       buildPhase = "true";
       installPhase = ''
         mkdir -p $out/bin
         cp $src/mytool $out/bin/
         chmod +x $out/bin/mytool
       '';
     };
   }

This will make ``mytool`` available as a package in ``pkgs``. You must then register it under your packages wherever needed by calling ``pkgs.mytool``.

----

**Tip:** Use this folder only for packages that cannot be obtained or easily overridden from upstream
