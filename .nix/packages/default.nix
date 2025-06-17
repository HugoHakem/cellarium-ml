# packages/default.nix
{ pkgs }:

let
  # List all .nix files in the directory (excluding this one)
  packageFiles = builtins.filter
    (file: file != "default.nix")
    (builtins.attrNames (builtins.readDir ./.));

  # Import each package file and merge their outputs into one attrset
  packageAttrs = builtins.foldl'
    (acc: file:
      acc // import (./. + "/${file}") { inherit pkgs; })
    {}
    packageFiles;
in
  packageAttrs
