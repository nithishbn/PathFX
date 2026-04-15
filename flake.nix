{
  description = "PathFX — drug-phenotype network enrichment analysis";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/b134951a4c9f3c995fd7be05f3243f8ecd65d798"; # nixos-24.05
  };

  outputs = {
    self,
    nixpkgs,
  }: let
    system = "x86_64-linux";
    pkgs = nixpkgs.legacyPackages.${system};

    fastcluster = pkgs.python311Packages.buildPythonPackage {
      pname = "fastcluster";
      version = "1.2.6";
      src = pkgs.fetchurl {
        url = "mirror://pypi/f/fastcluster/fastcluster-1.2.6.tar.gz";
        sha256 = "sha256-qriG76e2u6esEk9EmBU9BT5aCLgi0iVJJrcgbN9aiqY=";
      };
      buildInputs = [ pkgs.python311Packages.numpy ];
      doCheck = false;
    };

    myPython = pkgs.python311.withPackages (ps:
      with ps; [
        # Jupyter (optional, for interactive exploration)
        jupyter-core
        jupyterlab
        # PathFX core dependencies
        numpy
        pandas
        scipy
        matplotlib
        seaborn
        networkx
        fastcluster
        statsmodels
        openpyxl
        scikit-learn
      ]);

    # Text::OverlapFinder is bundled inside Text-Similarity — not in nixpkgs
    textSimilarity = pkgs.perlPackages.buildPerlPackage {
      pname = "Text-Similarity";
      version = "0.13";
      src = pkgs.fetchurl {
        url = "mirror://cpan/authors/id/T/TP/TPEDERSE/Text-Similarity-0.13.tar.gz";
        sha256 = "0v5d0gl5j82mk150nxn05j712wc7gxf4lpj7brh9xapa6iiqvhgr";
      };
      doCheck = false;
    };

    # UMLS::Interface — Perl interface to the UMLS MySQL database
    umlsInterface = pkgs.perlPackages.buildPerlPackage {
      pname = "UMLS-Interface";
      version = "1.51";
      src = pkgs.fetchurl {
        url = "mirror://cpan/authors/id/B/BT/BTMCINNES/UMLS-Interface-1.51.tar.gz";
        sha256 = "1nzl42l6arvjxk6hvrr3far0w9v9dyhdy3kw75iy6i2jqw05jcmm";
      };
      propagatedBuildInputs = with pkgs.perlPackages; [
        DigestSHA1
        DBI
        DBDmysql
      ];
      # Tests require a live UMLS MySQL instance
      doCheck = false;
    };

    # UMLS::Similarity — semantic similarity measures (provides umls-similarity.pl)
    umlsSimilarity = pkgs.perlPackages.buildPerlPackage {
      pname = "UMLS-Similarity";
      version = "1.47";
      src = pkgs.fetchurl {
        url = "mirror://cpan/authors/id/B/BT/BTMCINNES/UMLS-Similarity-1.47.tar.gz";
        sha256 = "0q7grj0qi50y7qpdp7fwyrh44hjqv4ly7fr2n43mdh0qgpmhjq8y";
      };
      propagatedBuildInputs = with pkgs.perlPackages; [
        umlsInterface
        textSimilarity
        LinguaStem
        TextNSP
      ];
      # Tests require a live UMLS MySQL instance
      doCheck = false;
    };

    myPerl = pkgs.perl.withPackages (_: [
      umlsInterface
      umlsSimilarity
    ]);

  in {
    devShells.x86_64-linux.default = pkgs.mkShell {
      buildInputs = [myPython myPerl];
      shellHook = ''
        echo "PathFX environment ready."
        echo "Run scripts from the scripts/ directory:"
        echo "  cd scripts && python phenotype_enrichment_pathway_Pfx050120.py -d <drug> -a <name>"
        echo ""
        echo "For UMLS clustering: configure rscs/itfc_confg.txt and ensure MySQL is running."
      '';
    };
  };
}
