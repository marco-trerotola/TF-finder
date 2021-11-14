#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# prototipi
sub Trace($);
sub fetchBond($);

# alcune variabili e array di uso generale
my $debug = 1;
my @aminoacid = qw/ALA ARG ASP ASN CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER TRP THR TYR VAL/;
my @nucAcid = qw/A T G C/;


my $bondFile = shift;
# scrive su file .outD e .outI
my ($outFile) =  split/\./, $bondFile;
my $outFileD = $outFile . '.outD';
my $outFileI = $outFile . '.outI';

my $lines = fetchBond($bondFile);	# contiene le linee del bond, rispulite da commenti
my @linee = @$lines;	# contiene le linee del bond, rispulite da commenti

Trace("\ncreo %hashCount...\n");
my %hashCount;
my %hashContatti;

Trace("\nInizio calcolo conteggio contatti");
foreach my $row (@$lines) {
    Trace("\n$row");
    my @cols = split(/\s+/, $row);
    if ( grep(/^$cols[0]$/, @aminoacid) and grep(/^$cols[4]$/, @nucAcid) ) {
	# like ASN A 31 ND2 A _ _ O2P _
	Trace("\n\tcontatto amminoacido-acidonucleico -> ($cols[0],$cols[1],$cols[2])++");
	my $chiave = $cols[0]."_".$cols[1]."_".$cols[2];
	Trace("\n\tchiave conteggio -> ($chiave)++");
	$hashCount{$chiave}++;
	# aggiungo il contatto
	if (($cols[7] =~ /O\dP/) or ($cols[7] eq 'P')) {
	    # like ASN A 31 ND2 A U 4 O1P 2.28 (fosfato)
	    Trace("\n\taggiungo a $chiave contatto OnP (fosfato, $cols[7])");
	    push @{$hashContatti{$chiave}}, 'OnP';
	} elsif ($cols[7] =~ /\*$/) {
	    # like SER A 65 O T R 7 C1* 3.09 (residuo)
	    Trace("\n\taggiungo a $chiave contatto Ribosio ($cols[7])");
	    push @{$hashContatti{$chiave}}, 'Ribosio';
	} else {
	    # like LYS A 32 NZ G S 13 N7 2.68
	    Trace("\n\taggiungo a $chiave contatto con acidonucleico ($cols[4])");
	    push @{$hashContatti{$chiave}}, $cols[4];
	}
    } elsif ( grep(/^$cols[0]$/, @aminoacid) and ($cols[4] eq 'HOH') ) {
	# like ARG C 336 NH1 HOH _ _ _
	Trace("\n\tPonte ($cols[0],$cols[4])");
	my $sublines = fetchBond($bondFile);	# contiene le linee del bond,
						# rispulite da commenti
	foreach my $subRow (@$sublines) {
	    my @subCols = split(/\s+/, $subRow);
	    if ( ( ($cols[5] =~ /^\d+$/) and ($subCols[0] eq 'HOH') and ($subCols[1] eq $cols[5]) and (!(grep(/^$subCols[3]$/, @aminoacid)))) or
		($cols[5] !~ /^\d+$/) and ($subCols[0] eq 'HOH') and ($subCols[1] eq $cols[5]) and ($subCols[2] eq $cols[6]) and (!(grep(/^$subCols[4]$/, @aminoacid)))) {
		# like	HIS A 328 ND1 HOH 1061 _ _ <--> HOH 1061 O A _ _ _ _
		# like HIS A 328 ND1 HOH B 1061 _ _ <--> HOH B 1061 O A _ _ _ _
		Trace("\n\t\t$subRow");
		Trace("\n\tcontatto amminoacido-HOH-acidonucleico -> ($cols[0],$cols[1],$cols[2])++");
		my $chiave = $cols[0].'_'.$cols[1].'_'.$cols[2];
		Trace("\n\tchiave conteggio -> ($chiave)++");
		$hashCount{$chiave}++;
		my $i;		# spiazzamento
		$i = $cols[5] =~ /^\d+$/ ? 0 : 1;
		# aggiungo il contatto
		if (($subCols[6+$i] =~ /O\dP/) or ($subCols[6+$i] eq 'P')) {
		    # like HOH ? 105 O C R 2 O2P 3.20
		    Trace("\n\taggiungo a $chiave contatto OnP (fosfato, $subCols[6+$i])");
		    push @{$hashContatti{$chiave}}, 'I_OnP';
		} elsif ($subCols[6+$i] =~ /\*$/) {
		    # like SER ? 65 O T R 7 C1* 3.09 (residuo)
		    Trace("\n\taggiungo a $chiave contatto Ribosio ($subCols[6+$i])");
		    push @{$hashContatti{$chiave}}, 'I_Ribosio';
		} else {
		    # like LYS ? 32 NZ G S 13 N7 2.68
		    Trace("\n\taggiungo a $chiave contatto con acidonucleico ($subCols[3+$i])");
		    push @{$hashContatti{$chiave}}, 'I_'.$subCols[3+$i];
		}
	    } else {
		Trace ".";	# debug con un puntino
	    }
	}
    } elsif (($cols[1] =~ /^\d+$/) and ($cols[0] eq 'HOH') and (grep(/^$cols[3]$/, @aminoacid))) {
	# like HOH 2 O LYS _ _ _ _
	Trace("\n\tPonte invertito ($cols[0],$cols[3])");
	my $sublines = fetchBond($bondFile);	# contiene le linee del bond,
						# rispulite da commenti
	foreach my $subRow (@$sublines) {
	    my @subCols = split(/\s+/, $subRow);
	    if ( (($subCols[0] eq 'HOH') and ($subCols[1] eq $cols[1]) and (!(grep(/^$subCols[3]$/, @aminoacid))) or
	    	 (($subCols[5] eq 'HOH') and ($subCols[6] eq $cols[1])) )) {
		# like HOH 2 O LYS A  275 _ _ <--> HOH 2 O A
		# like HOH 2 O LYS A  275 _ _ <--> _ A C 2 N6 HOH 2 _ _ _
		Trace("\n\t\t$subRow");
		Trace("\n\tcontatto inverso amminoacido-HOH-acidonucleico -> ($cols[3],$cols[4],$cols[5])++");
		my $chiave = $cols[3].'_'.$cols[4].'_'.$cols[5];
		Trace("\n\tchiave conteggio -> ($chiave)++");
		$hashCount{$chiave}++;
		# aggiungo il contatto
		if ($subCols[5] eq 'HOH') {
		    # like HOH ? 2 O LYS A  275 _ _ <--> _ A C 2 N6 HOH ? 2 _ _ _
		    if (($subCols[4] =~ /O\dP/) or ($subCols[4] eq 'P')) {
			# like _ A C 12 O2P HOH Z 42 O 3.13
			Trace("\n\taggiungo a $chiave contatto OnP (fosfato, $subCols[4])");
			push @{$hashContatti{$chiave}}, 'I_OnP';
		    } elsif ($subCols[4] =~ /\*$/) {
			# like _ A C 12 C1* HOH Z 42 O 3.13 (residuo)
			Trace("\n\taggiungo a $chiave contatto Ribosio ($subCols[4])");
			push @{$hashContatti{$chiave}}, 'I_Ribosio';
		    } else {
			# like _ A C 12 N7 HOH Z 42 O 3.13
			Trace("\n\taggiungo a $chiave contatto con acidonucleico ($subCols[1])");
			push @{$hashContatti{$chiave}}, 'I_'.$subCols[1];
		    }
		} else {
		    # like HOH ? 2 O LYS A  275 _ _ <--> HOH ? 2 O A
		    my $i;		# spiazzamento
		    $i = $cols[1] =~ /^\d+$/ ? 0 : 1;
		    if (($subCols[6+$i] =~ /O\dP/) or ($subCols[6+$i] eq 'P')) {
			# like HOH ? 105 O C R 2 O2P 3.20
			Trace("\n\taggiungo a $chiave contatto OnP (fosfato, $subCols[6+$i])");
			push @{$hashContatti{$chiave}}, 'I_OnP';
		    } elsif ($subCols[6+$i] =~ /\*$/) {
			# like HOH ? 105 O C R 2 C1* 3.20
			Trace("\n\taggiungo a $chiave contatto Ribosio ($subCols[6+$i])");
			push @{$hashContatti{$chiave}}, 'I_Ribosio';
		    } else {
			# like HOH ? 105 O C R 2 N7 3.20
			Trace("\n\taggiungo a $chiave contatto con acidonucleico ($subCols[3+$i])");
			push @{$hashContatti{$chiave}}, 'I_'.$subCols[3+$i];
		    }
		}
	    } else {
		Trace ".";	# debug con un puntino
	    }
	}
    } elsif (($cols[2] =~ /^\d+$/) and ($cols[0] eq 'HOH') and (grep(/^$cols[4]$/, @aminoacid))) {
	# like HOH Z 2 O LYS _ _ _ _
	Trace("\n\tPonte invertito ($cols[0],$cols[4])");
	my $sublines = fetchBond($bondFile);	# contiene le linee del bond,
						# ripulite da commenti
	foreach my $subRow (@$sublines) {
	    my @subCols = split(/\s+/, $subRow);
	    if ( ( ($cols[5] =~ /^\d+$/) and ($subCols[0] eq 'HOH') and ($subCols[1] eq $cols[5]) and (!(grep(/^$subCols[3]$/, @aminoacid)))) or
		($subCols[0] eq 'HOH') and ($subCols[1] eq $cols[5]) and ($subCols[2] eq $cols[6]) and (!(grep(/^$subCols[4]$/, @aminoacid)))) {
		# like HIS A 328 ND1 HOH 1061 _ _ <--> HOH 1061 O A _ _ _ _
		# like HIS A 328 ND1 HOH B 1061 _ _ <--> HOH B 1061 O A _ _ _ _
		Trace("\n\t\t$subRow");
		Trace("\n\tcontatto amminoacido-HOH-acidonucleico -> ($cols[0],$cols[1],$cols[2])++");
		my $chiave = $cols[0].'_'.$cols[1].'_'.$cols[2];
		Trace("\n\tchiave conteggio -> ($chiave)++");
		$hashCount{$chiave}++;
		# like HIS A 328 ND1 HOH ? 1061 _ _ <--> HOH ? 2 O A
		my $i;		# spiazzamento
		$i = $cols[5] =~ /^\d+$/ ? 0 : 1;
		if (($subCols[6+$i] =~ /O\dP/) or ($subCols[6+$i] eq 'P')) {
		    # like HOH ? 105 O C R 2 O2P 3.20
		    Trace("\n\taggiungo a $chiave contatto OnP (fosfato, $subCols[6+$i])");
		    push @{$hashContatti{$chiave}}, 'I_OnP';
		} elsif ($subCols[6+$i] =~ /\*$/) {
		    # like HOH ? 105 O C R 2 C1* 3.20
		    Trace("\n\taggiungo a $chiave contatto Ribosio ($subCols[6+$i])");
		    push @{$hashContatti{$chiave}}, 'I_Ribosio';
		} else {
		    # like HOH ? 105 O C R 2 N7 3.20
		    Trace("\n\taggiungo a $chiave contatto con acidonucleico ($subCols[3+$i])");
		    push @{$hashContatti{$chiave}}, 'I_'.$subCols[3+$i];
		}
	    } else {
		Trace ".";	# debug con un puntino
	    }
	}
    } elsif ( grep(/^$cols[5]$/, @aminoacid) and grep(/^$cols[1]$/, @nucAcid) ) {
	# like _ C C 5 N4 GLU A 9 OE1 3.04
	Trace("\n\tcontatto amminoacido-acidonucleico invertito -> ($cols[5],$cols[6],$cols[7])++");
	my $chiave = $cols[5]."_".$cols[6]."_".$cols[7];
	Trace("\n\tchiave conteggio -> ($chiave)++");
	$hashCount{$chiave}++;
	# aggiungo il contatto
	if (($cols[4] =~ /O\dP/) or ($cols[4] eq 'P')) {
	    # like _ A R 4 O2P GLN A 27 OE1 3.00
	    Trace("\n\taggiungo a $chiave contatto OnP (fosfato, $cols[4])");
	    push @{$hashContatti{$chiave}}, 'OnP';
	} elsif ($cols[4] =~ /\*$/) {
	    # like _ A R 4 O2* GLN A 27 OE1 3.00
	    Trace("\n\taggiungo a $chiave contatto Ribosio ($cols[4])");
	    push @{$hashContatti{$chiave}}, 'Ribosio';
	} else {
	    # like _ A R 4 N6 GLN A 27 OE1 3.00
	    Trace("\n\taggiungo a $chiave contatto con acidonucleico ($cols[1])");
	    push @{$hashContatti{$chiave}}, $cols[1];
	}
    } elsif ( grep(/^$cols[1]$/, @nucAcid) and ($cols[5] eq 'HOH') ) {
	# like _ C C 7 O3* HOH 72 O 2.60
	Trace("\n\tprobabile ponte amminoacido-acidonucleico invertito:");
	Trace("\n\ti ponti li attraverso dall'altra estremita'");
    } elsif ( (($cols[1] =~ /^\d+$/) and (grep(/^$cols[3]$/, @nucAcid)) and ($cols[0] eq 'HOH') ) xor
    		(($cols[2] =~ /^\d+$/) and (grep(/^$cols[4]$/, @nucAcid)) and ($cols[0] eq 'HOH')) ) {
	# like HOH 65 O G C 17 O1P 2.34
	# like HOH Z 7 O A _ _ _ _
	Trace("\n\tprobabile ponte contatto amminoacido-acidonucleico:");
	Trace("\n\ti ponti li attraverso dall'altra estremita'");
    } else {
	Trace("\n\n\tAttenzione: _NON_ contemplato!!\n");
    }
}

my $tot_contatti;

print "\n", Dumper(\%hashCount);
$tot_contatti = 0;
foreach $_ (values %hashCount) {
    $tot_contatti += $_;
}
Trace("\nTotale contatti: $tot_contatti");
Trace("\nFine calcolo conteggio contatti");

print "\n",Dumper(\%hashContatti);
$tot_contatti = 0;
foreach my $array (values %hashContatti) {
    foreach $_ (@$array) {
	$tot_contatti += 1;
    }
}
Trace("\nTotale contatti: $tot_contatti");
Trace("\nverificare che questa voce ed il conteggio precedente siano uguali");

Trace("\ncreo l'hashSintesiCount (hash frazionario)...");
my %hashSintesiCount = ();
foreach $_ (keys %hashCount) {
    my $chiave = $_;
    $hashSintesiCount{$chiave} += (int(1/$hashCount{$chiave}*100)/100);
}
Trace(" fine sintesi frazionaria\n");
print "\n", Dumper(\%hashSintesiCount);

Trace("\nAzzero l'hash per la sintesi finale... ");
my %hash = ();
foreach my $amnisD (@aminoacid)  {
    my $amnisI = $amnisD . '_I';
    foreach $_ (qw /OnP Ribosio A T G C/) {
	$hash{$amnisD}{$_} = 0;
	$hash{$amnisI}{$_} = 0;
    }
}
Trace('fine azzermento hash di sintesi finale');
#print "\n", Dumper(\%hash);

Trace("\naggiornamento hash finale...");
foreach my $key_row (keys %hashContatti) {
    my ($key) = split(/_/,$key_row);		# like (GLU) = split /_/, GLU_B_332
    foreach $_ (@{$hashContatti{$key_row}}) {
        if ($_ =~ /^I_/) {
	    $_ =~ s/^I_(.*)$/$1/;
	    $hash{$key.'_I'}{$_} += $hashSintesiCount{$key_row};
	} else {
	    $hash{$key}{$_} += $hashSintesiCount{$key_row};
	}
    }
}
Trace("fine aggiornamento hash finale");
print "\n", Dumper(\%hash);

Trace("\n\nSintesi di $bondFile in $outFileD e $outFileI");
open(wD, ">$outFileD");
open(wI, ">$outFileI");
Trace("\nSintesi diretta...");
#printf wD ("%6s %5s %5s %5s %5s %5s %5s\n", '','A','T','G','C','Fosfato','Ribosio');
printf wD "\tA\tT\tG\tC\tFosfato\tRibosio";
#print wD '-' x45, "\n";
foreach my $aa (sort @aminoacid) {
#    printf wD ("%6s|%5s|%5s|%5s|%5s|%5s|%5s\n", $aa, $hash{$aa}{'A'}, $hash{$aa}{'T'}, $hash{$aa}{'G'}, $hash{$aa}{'C'}, $hash{$aa}{'OnP'}, $hash{$aa}{'Ribosio'});
    print wD "\n$aa\t$hash{$aa}{'A'}\t$hash{$aa}{'T'}\t$hash{$aa}{'G'}\t$hash{$aa}{'C'}\t$hash{$aa}{'OnP'}\t$hash{$aa}{'Ribosio'}";
}
Trace(" ok\n");
Trace("Sintesi mediata da HOH...");
#printf wI ("%6s|%5s|%5s|%5s|%5s|%5s|%5s\n", '','A','T','G','C','Fosfato','Ribosio');
printf wI "\tA\tT\tG\tC\tFosfato\tRibosio";
#print wI '-' x45, "\n";
foreach my $aa (sort @aminoacid) {
    $aa = $aa . "_I";
#    printf wI ("%6s|%5s|%5s|%5s|%5s|%5s|%5s\n", $aa, $hash{$aa}{'A'}, $hash{$aa}{'T'}, $hash{$aa}{'G'}, $hash{$aa}{'C'}, $hash{$aa}{'OnP'}, $hash{$aa}{'Ribosio'});
    print wI "\n$aa\t$hash{$aa}{'A'}\t$hash{$aa}{'T'}\t$hash{$aa}{'G'}\t$hash{$aa}{'C'}\t$hash{$aa}{'OnP'}\t$hash{$aa}{'Ribosio'}";
}
Trace(" ok\n\n\n");
close (wD);
close (wI);



sub fetchBond($) {
    my $bondFile = shift;
    die("bond $bondFile non trovato") unless -e "$bondFile";
    open (RH, "<$bondFile");
    Trace("\nFetch bond...");
    my @lines = ();
    while (<RH>) {
	my $row = $_;
	chop $row;
	my @colonne = split(/\s+/, $row);			# non la considero
	unless ($#colonne >=6) {				# se la linea 
	    # Trace("$row -> salto < 6 colonne\n");	# ha meno di
	    next					# 6 colonne
	}
	push @lines, $row;
    }
    close RH;
    Trace("fine fetch bond...");
    return \@lines;
}

sub Trace($) {
    ($_) = @_;
    print $_ if $debug == 1;
}

