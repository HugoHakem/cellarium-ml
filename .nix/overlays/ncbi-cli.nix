self: super: 
let
  inherit (super.stdenv.hostPlatform) system;
  platform = if system == "x86_64-linux" then "linux-amd64"
             else if system == "aarch64-linux" then "linux-arm64"
             else throw "Unsupported system: ${system}";
  baseUrl = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/${platform}";
in {
  ncbi-datasets-cli = super.stdenv.mkDerivation {
    pname = "ncbi-datasets-cli";
    version = "2.0.0";
    src = super.fetchurl {
      url = "${baseUrl}/datasets";
      sha256 = "sha256-+RRXk5t8O98jC0jcfbY+eQWbz1ulZNa4RCNHXS70804=";
    };
    phases = [ "installPhase" ];
    installPhase = ''
      mkdir -p $out/bin
      cp $src $out/bin/datasets
      chmod +x $out/bin/datasets
    '';
  };
  ncbi-dataformat-cli = super.stdenv.mkDerivation {
    pname = "ncbi-dataformat-cli";
    version = "2.0.0";
    src = super.fetchurl {
        url = "${baseUrl}/dataformat";
        sha256 = "sha256-ZGOQss8e9jXCZxUXpkDntu0rETroFqkwBBnAWKVddxE";
    };
    phases = [ "installPhase" ];
    installPhase = ''
        mkdir -p $out/bin
        cp $src $out/bin/dataformat
        chmod +x $out/bin/dataformat
    '';
    };
}