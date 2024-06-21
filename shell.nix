{ pkgs ? import <nixpkgs> { } }:
with pkgs;
let
  simquantgen = pkgs.rPackages.buildRPackage {
    name = "simquantgen";
    src = pkgs.fetchFromGitHub {
      owner = "jeffersonfparil";
      repo = "simquantgen";
      rev = "fed3612a22afa70e3bea758addbe0da421fad7ce";
      sha256 = "sha256-ZGGKV7eR6u7ll4ulejyS0hoji20s4tIjqiyncalZ3Og=";
    };
    propagatedBuildInputs = with rPackages; [
        MASS
        doParallel
        foreach
        txtplot
    ];
  }; ### Use lib.fakeSha256 on "sha256" to get the correct code from the error messages
  my-r-pkgs = rWrapper.override {
    packages = with rPackages; [
      remotes
      devtools
      testthat
      vcfR
      txtplot
      glmnet
      BGLR
      sommer
      MASS
      parallel
      doParallel
      foreach
      shiny
      shinyWidgets
      shinyFiles
      shinycssloaders
      plotly
      bslib
      DBI
      RSQLite
      simquantgen
    ];
  };
in mkShell {
  buildInputs = with pkgs; [git glibcLocales openssl which openssh curl wget ];
  inputsFrom = [ my-r-pkgs ];
  shellHook = ''
    mkdir -p "$(pwd)/_libs"
    export R_LIBS_USER="$(pwd)/_libs"
  '';
  GIT_SSL_CAINFO = "${cacert}/etc/ssl/certs/ca-bundle.crt";
#  LOCALE_ARCHIVE = stdenv.lib.optionalString stdenv.isLinux
#    "${glibcLocales}/lib/locale/locale-archive";
}

### nix-shell --run bash --pure