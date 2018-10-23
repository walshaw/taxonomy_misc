#!/usr/bin/perl

use strict;
use warnings;

use Carp;

use Getopt::Long;
use List::Util qw( first );

my $delimiter        = q{|};
my $delimiter_quoted = quotemeta $delimiter;

my $taxonomy_names_file;
my $nodes_list_file;
my $input_fasta_file;
my $case_sensitive; # N.B. *only* applies to the comparison of taxon strings (read
                    # from FASTA file) and contents of $taxonomy_names_file; in all
                    # other contexts, *the original case is always used*
my $resolve_higher;
my $report_frequency = 20000;

my $usage = qq{Usage:\n\n$0  [ options ]

Options:
        -namesfile	 string		TAXONOMY_NAMES_FILE
	-nodeslistfile	 string		NODES_LIST_FILE
	-fastafile	 string		FASTA_FILE
        -casesensitive   		(boolean; default is false)
        -resolvehigher                  (boolean; default is false)
        -reportfrequency integer	REPORT_FREQUENCY (default $report_frequency)

-namesfile , -nodeslistfile and -fastafile are mandatory arguments; if the switches
are omitted, then the first, second and third unswitched arguments are interpreted
respectively as TAXONOMY_NAMES_FILE, NODES_LIST_FILE, FASTA_FILE.

TAXONOMY_NAMES_FILE should be NCBI taxonomy's names.dmp, or something very similar.

NODES_LIST_FILE is a file containing the IDs of all the nodes required, one node ID
per line (only the first column is read). Such a file can be obtained by extracting
the second field from a file output by ncbi_taxonomy_clade.pl (example lines from
such a file, where the object was to identify all nodes of the proteobacteria clade,
whose ancestral node has ID 1224, is as follows):

node        1   IN CLADE:NO     terminal node        1; lineage: 1
node        2   IN CLADE:NO     terminal node        1; lineage: 2, 131567, 1
node        6   IN CLADE:YES    terminal node     1224; lineage: 6, 335928, 356, 28211, 1224
node        7   IN CLADE:YES    terminal node     1224; lineage: 7, 6
node        9   IN CLADE:YES    terminal node     1224; lineage: 9, 32199, 543, 91347, 1236, 1224

FASTA_FILE must have headers which end in the bracket-enclosed taxonomic strings, i.e.
the NCBI format, e.g.:

>gi|418280395|ref|ZP_12893276.1| hypothetical protein SA21178_2416 [Staphylococcus aureus subsp. aureus 21178]
>gi|260642595|ref|ZP_05859457.1| conserved hypothetical protein [Bacteroides finegoldii DSM 17565]
>gi|336319454|ref|YP_004599422.1| hypothetical protein Celgi_0332 [[Cellvibrio] gilvus ATCC 13127]
>gi|526218031|ref|YP_008263104.1| purine biosynthesis protein purH [Salmonella enterica subsp. enterica serovar 4,[5],12:i:- str. 08-1736]
>gi|32477681|ref|NP_870675.1| gluconolactonase [precursor] [Rhodopirellula baltica SH 1]
>gi|218132225|ref|ZP_03461029.1| hypothetical protein BACPEC_00082 [[Bacteroides] pectinophilus ATCC 43243]

- note that pairs of brackets can be nested, and also the '[precursor]' construct. The vast
majority of the NCBI data has the simple format as in the first two example lines above.

By default, when comparing a taxon name to those in the TAXONOMY_NAMES_FILE, a
case-insensitive match is done; this can cause problems e.g. the case of 'Vibrio sp. Ex25'
which appears in the names.dmp file only as 'Vibrio sp. EX25'. If -casesensitive is
specified, then the original case is compared in all cases.

If there is a complete failure to find a taxon name in the list, then -resolvehigher will
invoke an attempt (in those cases only) to use progressively higher taxonomic groups. E.g.
currently, this taxon string appears in NC_016610.faa: 'Tannerella forsythia ATCC 43037'
- but is currently absent from names.dmp, which does however contain an entry for
'Tannerella forsythia'. Basically, if -resolvehigher is specified, then where a failed
entry contains N subfields (N=4 in this example), progressively fewer will be used for
the next test (i.e. first 'Tannerella forsythia ATCC', then 'Tannerella forsythia',
then (if necessary - not in this example), 'Tannerella' would be tested) until the key
is found or no subfields remain (in which case, it's fatal).


};

croak qq{$usage} if !(GetOptions(
        "namesfile=s"       => \$taxonomy_names_file,
        "nodeslistfile=s"   => \$nodes_list_file,
        "fastafile=s"       => \$input_fasta_file, # currently, output FASTA is sent to STDOUT
        "casesensitive"     => \$case_sensitive,
        "resolvehigher"     => \$resolve_higher,
        "reportfrequency=i" => \$report_frequency,
    )
);

# these are mandatory arguments
$taxonomy_names_file ||= shift;
$nodes_list_file     ||= shift;
$input_fasta_file    ||= shift;

croak $usage if !(defined $taxonomy_names_file && defined $nodes_list_file && defined $input_fasta_file);

my %name_is_required; # keys are taxon names; values, if defined, mean:
                      #    non-zero means, it's known that this taxon is required
                      #    zero means, it's known that this taxon is NOT required
                      # - whereas if undef, then the true/false value must be
                      # resolved and set

print STDERR qq{reading node list file $nodes_list_file\n\n};
my $required_nodes_lref  = read_node_list($nodes_list_file);

print STDERR qq{reading taxonomy names file $taxonomy_names_file\n\n};
my ($node_id_of_href)    = read_names($taxonomy_names_file, $case_sensitive);

#my @required_nodes = @{required_nodes_lref};
#my @node_id_of     = @{$node_id_of_href};

print STDERR qq{reading FASTA file $input_fasta_file\n\n};

process_fasta( fasta             => $input_fasta_file   ,
               required_node_ids => $required_nodes_lref,
               node_ids          => $node_id_of_href    ,
               required_names    => \%name_is_required  ,
               case_sensitive    => $case_sensitive     ,
               resolve_higher    => $resolve_higher     ,
               report_frequency  => $report_frequency   ,
             );

exit 0;

sub read_node_list {

    my $nodes_list_file = shift or croak qq{supply nodes list file name};

    my @required_nodes;

    open my $fh, '<', $nodes_list_file or croak qq{couldn't read nodes list file $nodes_list_file};

    while (defined (my $line = <$fh>)) {
        chomp $line;
        $line =~ s{ \A \s* (\d+) \s* (?: \S .* )? \z }{$1}xms; # get rid of extra whitespace and any extra columns
        push @required_nodes, $line;
    }

    close $fh or croak qq{couldn't close nodes list file $nodes_list_file};

    return \@required_nodes;
}


sub read_names {

    my $taxon_names_file = shift or croak qq{supply taxonomy names/node IDs file name};
    my $case_sensitive   = shift;

    my %node_id_of; # keys are taxon names - as in second field below ($node_name); values are node IDs
    my @node_names; # indices represent node IDs; many won't be used; values are the names (second field
                    # - $node_name - see below) of rows where the fourth field ($class) is
                    # 'scientific name' - other rows are ignored (otherwise there wouldn't be a 1:1
                    # ID:name mapping)

    open my $fh, '<', $taxon_names_file or croak qq{couldn't read taxonomy names/node IDs file $taxon_names_file};

    while (defined (my $line = <$fh>)) {
        chomp $line;
        my @fields = split m{ \t $delimiter_quoted \t }xms, $line; # note the explicit tabs; using \s* would
                                                                   # cause a minor problem with null field values

        # see ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt
        # - the four fields' official names are tax_id, name_txt, 'unique name' (sic), 'name class' (sic)

        my ($node_id, $node_name, $unique_name, $class) = @fields;

        my $comparison_node_name = $node_name;
        if (!$case_sensitive) { $comparison_node_name = lc $node_name; }

        $node_id_of{$comparison_node_name} = $node_id;

        if ($class eq 'scientific name') {
            $node_names[$node_id] = $node_name; # note that this uses the original name regardless of case-sensitivity
        }

    }

    close $fh or croak qq{couldn't close taxonomy names/node IDs file $taxon_names_file};

    return \%node_id_of, \@node_names;
}


sub process_fasta {

    my %arg = @_;

    my ($input_fasta_file      , $required_nodes_lref, $node_id_of_href,
        $name_is_required_href , $case_sensitive     , $resolve_higher, $report_frequency)
      = @arg{'fasta','required_node_ids','node_ids','required_names','case_sensitive', 'resolve_higher', 'report_frequency'};

    my @current_fasta_record;
    ###my @next_fasta_record;

    croak qq{supply FASTA file name} if (!$input_fasta_file);

    # decided not to use BioPerl routines for this very simple FASTA parsing;
    # want to make it as fast as possible, so keep overheads to a minimum

    open my $fh, '<', $input_fasta_file or croak qq{couldn't read input FASTA file $input_fasta_file};

    my $n_records = 0; # total, irrespective of whether used/output or not
    my $n_lines   = 0;

    my $n_positive_records = 0;

    while (defined (my $line = <$fh>)) {

        $n_lines++;

        if ($line =~ m{ \A [>] }xms) {

            if (@current_fasta_record) {
                my ($seq_id, $description, $taxon) = process_fasta_record(\@current_fasta_record);

                my $result = test_taxon(
                    taxon             => $taxon                 ,
                    required_node_ids => $required_nodes_lref   ,
                    node_ids          => $node_id_of_href       ,
                    required_names    => $name_is_required_href ,
                    case_sensitive    => $case_sensitive        ,
                    resolve_higher    => $resolve_higher        ,
                );

                if ($result) { output_fasta_record(\@current_fasta_record); $n_positive_records++; }
            }

            $n_records++;
            if ($report_frequency && !($n_records % $report_frequency)) {
                print STDERR qq{$n_records records processed; $n_positive_records are positive; $n_lines lines read\n};
            }
            @current_fasta_record = ($line);
        }
        else {
            push @current_fasta_record, $line;
        }

    }

    if (@current_fasta_record) {

        my ($seq_id, $description, $taxon) = process_fasta_record(\@current_fasta_record);

        my $result = test_taxon(
            taxon             => $taxon                 ,
            required_node_ids => $required_nodes_lref   ,
            node_ids          => $node_id_of_href       ,
            required_names    => $name_is_required_href ,
            case_sensitive    => $case_sensitive        ,
            resolve_higher    => $resolve_higher        ,
            );

        if ($result) { output_fasta_record(\@current_fasta_record); $n_positive_records++; }

        $n_records++;
    }


    print STDERR qq{total of $n_lines lines and $n_records records read ($n_positive_records are positive) from FASTA file $input_fasta_file\n\n};

    close $fh or croak qq{couldn't close input FASTA file $input_fasta_file};


}

sub test_taxon {

    my %arg = @_;

    my ($taxon,  $required_nodes_lref,  $node_id_of_href,  $name_is_required_href, $case_sensitive, $resolve_higher)
      = @arg{'taxon','required_node_ids','node_ids','required_names','case_sensitive','resolve_higher'};

    # N.B. $case_sensitive *only* applies to the comparison of taxon strings (read
    # from FASTA file) and contents of $taxonomy_names_file; in all
    # other contexts, *the original case is always used*

    my $result = 0;

    # If $name_is_required_href itself is not defined, then this FASTA record
    # automatically tests 'positive'; i.e. no hashref means 'use everything'

    return 1 if (!defined $name_is_required_href);

    # If $name_is_required_href itself is defined, then the current taxon name
    # may or may not exist as a key; if it does, then the corresponding value
    # is the answer

    return $name_is_required_href->{$taxon} if (defined $name_is_required_href->{$taxon});

    # If $name_is_required_href->{$taxon} is undefined, then attempt to get
    # the right answer - which relies on $required_nodes_lref and $node_id_of_href
    # both being defined.

    # So if either  $required_nodes_lref or $node_id_of_href are NOT defined,
    # then again default to 'positive':

    return 1 if !(defined $required_nodes_lref && defined $node_id_of_href);

    # Ok - so here we have all we need; presumably. Problem is, if $node_id_of_href
    # IS defined by has no key corresponding to the current taxon, then it's a
    # sign that something is very wrong; so bail out.

    # There are some gotchas here. E.g. this file:

    # ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Vibrio_Ex25_uid41601/NC_013456.faa

    # Has FASTA headers where the taxon string is unambiguously 'Vibrio sp. Ex25';
    # yet there is no exact name matching this in the NCBI's current names.dmp; this
    # file:
    # ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Vibrio_Ex25_uid41601/NC_013456.rpt
    #
    # confirms the taxon name:
    #
    # Taxname: Vibrio sp. Ex25
    #
    # but checking the contents of the names.dmp file:
    # $ grep 'Vibrio' names.dmp  | grep -i Ex25
=pod
150340	|	Vibrio sp. EX25	|		|	includes	|
150340	|	Vibrio sp. Ex25_939	|		|	misspelling	|
$ grep '^150340   ' names.dmp # that's a tab after the number
150340	|	"Vibrio antiquarius" Hasan et al. 2015	|		|	authority	|
150340	|	Vibrio antiquarius	|		|	scientific name	|
150340	|	Vibrio sp. BV25Ex	|		|	misspelling	|
150340	|	Vibrio sp. EX25	|		|	includes	|
150340	|	Vibrio sp. Ex25_939	|		|	misspelling	|
150340	|	strain EX25	|		|	type material	|
=cut
    # - so for some reason the FASTA file cites the 'includes' field - and with
    # inconsistent case at that - rather than the scientific name 'Vibrio antiquarius'

    # Another reminder that the only time that case-(non-)sensitivity is an issue
    # is when setting/getting the value corresponding to the key in the
    # %{$node_id_of_href} hash

    ###my @taxon_subfields;
    ###my $taxon_truncated;

    if (!(defined $node_id_of_href->{$taxon} || defined $node_id_of_href->{lc $taxon})) {
        my $warning = qq{no key '$taxon' (or '}.lc($taxon).qq{') in hash of node names => IDs};
        croak qq{FATAL: $warning} if (!$resolve_higher || ( ($taxon !~ m{ \s }xms) && ($taxon !~ m{ [\[\]] }xms)) );
        print STDERR qq{WARNING: $warning\n};
#        $taxon_truncated =~ s{ \s+ \S+ \z }{}xms;
    }
    my $node_id;

    GET_NODE_ID: {
      if ($case_sensitive) {

        $node_id = $node_id_of_href->{$taxon};

        if (!defined $node_id) {
            my $warning = qq{no key '$taxon' in hash of node names => IDs};

            croak qq{FATAL: $warning} if (!$resolve_higher || ( ($taxon !~ m{ \s }xms) && ($taxon !~ m{ [\[\]] }xms) ) );

            # Another annoying exception: the taxon string:
            # appears in the NCBI FASTA headers:
            # '[Eubacterium] cylindroides T2-87'
            # but there is no such string in the current names.dmp, but only:
            # 'Eubacterium cylindroides T2-87'
            #
            # Note that this contrasts with many other similar format strings, e.g.:
            # '[Cellvibrio] gilvus ATCC 13127' which does appear literally
            # in names.dmp

            # So the first step after a failure to find the key, is to try a 'bracketless'
            # version

            (my $bracketless_taxon = $taxon) =~ s{ [\[\]] }{}gxms;
            print STDERR qq{\ttrying taxon string '$bracketless_taxon'\n};
            $node_id = $node_id_of_href->{$bracketless_taxon};
            if (defined $node_id) {
                print STDERR qq{WARNING: successfully resolved '$bracketless_taxon' to node $node_id; using as proxy for '$taxon'\n\n};
                last GET_NODE_ID;
            }

            # An example of a troublesome taxon string is 'Tannerella forsythia ATCC 43037'
            # - see more in the usage string above.

            print qq{WARNING: $warning\n};
            my $taxon_truncated = $taxon;
            while (!defined($node_id) && ($taxon_truncated =~ s{ \s+ \S+ \z }{}xms)) {
                print STDERR qq{\ttrying taxon string '$taxon_truncated'\n};
                $node_id = $node_id_of_href->{$taxon_truncated};
            }
            croak qq{coudln't resolve '$taxon'} if (!defined $node_id);
            print STDERR qq{\tresolved '$taxon_truncated' to node $node_id; using as proxy for '$taxon'\n\n};
=pod
            # recursive call; eventually, either success (the value of the
            # ancestral ) will be returned,
            # or one of the nested calls will result in a fatal error
            $result = test_taxon(
                       taxon             => $taxon_truncated       ,
                       required_node_ids => $required_nodes_lref   ,
                       node_ids          => $node_id_of_href       ,
                       required_names    => $name_is_required_href ,
                       case_sensitive    => $case_sensitive        ,
                       resolve_higher    => $resolve_higher        ,
                      );
=cut
        }
      }
      else {

        $node_id = $node_id_of_href->{lc $taxon};

        if (!defined $node_id) {
            my $warning = qq{$taxon: no key '}.lc($taxon).qq{' in hash of node names => IDs};
            croak qq{FATAL: $warning} if (!$resolve_higher || ( ($taxon !~ m{ \s }xms) && ($taxon !~ m{ [\[\]] }xms) ) );

            # First step after a failure to find the key, is to try a 'bracketless'
            # version

            (my $bracketless_taxon = lc $taxon) =~ s{ [\[\]] }{}gxms;
            print STDERR qq{\ttrying taxon string '$bracketless_taxon'\n};
            $node_id = $node_id_of_href->{$bracketless_taxon};
            if (defined $node_id) {
                print STDERR qq{WARNING: successfully resolved '$bracketless_taxon' to node $node_id; using as proxy for '$taxon'\n\n};
                last GET_NODE_ID;
            }

            print qq{WARNING: $warning\n};
            my $taxon_truncated = lc $taxon;
            while (!defined($node_id) && ($taxon_truncated =~ s{ \s+ \S+ \z }{}xms)) {
                print STDERR qq{\ttrying taxon string '$taxon_truncated'\n};
                $node_id = $node_id_of_href->{$taxon_truncated};
            }
            croak qq{coudln't resolve '$taxon'} if (!defined $node_id);
            print STDERR qq{\tresolved '$taxon_truncated' to node $node_id; using as proxy for '$taxon'\n\n};
        }
        # Note that the hash is such that there either the original-case version or the lowercase
        # version will be (is expected to be) present as a key, but not both
      }
    } # end GET_NODE_ID block

    if (defined first { $_ == $node_id } @{$required_nodes_lref}) {
###        $result = $node_id; # relies on there being no node 0
        $result = 1;
    }

    # Now set the result for this taxon name, so that this check doesn't have
    # to be repeated

    $name_is_required_href->{$taxon} = $result;

    return $result;

}

sub process_fasta_record {

    # decided not to use BioPerl routines for this very simple FASTA parsing

    my $fasta_record_lref = shift;

    my @fasta_record = @{$fasta_record_lref};

    my $header = $fasta_record[0];

    my ($seq_id, $description, $taxon) = parse_fasta_header($header);

    ###print qq{$seq_id\n\t'$header'\n'$description'\n\t'$taxon'\n\n};

    return ($seq_id, $description, $taxon);
}

sub parse_fasta_header {

    my $header = shift;

=pod
    example headers:

>gi|260642595|ref|ZP_05859457.1| conserved hypothetical protein [Bacteroides finegoldii DSM 17565]
>gi|336319454|ref|YP_004599422.1| hypothetical protein Celgi_0332 [[Cellvibrio] gilvus ATCC 13127]
>gi|526218031|ref|YP_008263104.1| purine biosynthesis protein purH [Salmonella enterica subsp. enterica serovar 4,[5],12:i:- str. 08-1736]
>gi|32477681|ref|NP_870675.1| gluconolactonase [precursor] [Rhodopirellula baltica SH 1]
>gi|260599358|ref|YP_003211929.1| [Citrate [pro-3S]-lyase] ligase [Cronobacter turicensis z3032]
>gi|386308790|ref|YP_006004846.1| gpe+E' [Enterobacteria phage P2] gb [AAD03292.1]
>gi|479170547|ref|YP_007798872.1| hypothetical protein [[Eubacterium] cylindroides T2-87]

=cut

    # Regexp is based closely on one I used in the parse_hmmer_domtable.pl script; this fails
    # due to the 5th example above:

#    croak qq{unexpected header:\n$header\n} if $header !~ m{ \A [>] (\S+) \s+ ( [^\[]+ [^\[\s] (?: \s* \[ precursor \] )? ) \s* \[ ( .+ ) \] \s* \z }xms;

    # the 6th example above is particularly hideous; it genuinely appears to be a one of a kind;
    # it is literally the only case of '...] gb [...]' in any of the complete, draft or plasmid
    # sets. It's the only case I'm aware of where the final '[...]' doesn't contain the taxon
    # string.

    my ($desc_remainder, $taxon_string);

    (my $description = $header) =~ s{ \A [>] (\S+) \s + (\S .* \S) \s* \n? \z }{$2}xms;
    my $seq_id = $1;

    croak qq{unexpected header:\n$header\n} if (!$seq_id);

    # deal with the horrible '6th example' above and any future wayward records in the same vein

    if ($description =~ s{ ( [\]] ) \s+ gb \s+ [\[] [A-Z0-9\.]+ [\]] \z }{$1}xms) {
        print STDERR qq{WARNING: removed spurious trailing text from this header:\n$header\n\n};
    }

    # strategy is to capture the easy ones first; saves time, because most of the lines will be like
    # this:

    if ( $description =~ m{ \A ([^\[\]]+) \s+ \[ ( [^\[\]]+ ) \] \s* \z }xms) {
        ($desc_remainder, $taxon_string) = ($1, $2);
    }
    else {
        # otherwise, strategy is to find the final ']' in the header, then simply identify the
        # corresponding '[' character; everything between the two (which *can in principle
        # include **any number** of nested pairs of '[' .. ']' *) is the taxon string;
        # everything before that corresponding '[' character is not (irrespective of whether
        # it contains any '[' .. ']' pairs, as they will be things like,
        # '[Citrate [pro-3S]-lyase]', '[precursor]' etc; see examples above.
        # Because of the unknown numbers of nested '[ .. ]' in the taxon string, a regexp
        # could be a bit hairy, so it's done a different way, simply working from the back:


        # (first do a check that the header ends in ']' and ensure this is the final character)
        croak qq{unexpected header:\n$header\n} if ($description !~ s{ ( [ \] ] ) \s* \n? \z }{$1}xms);

        my $index = length($description) - 2; # position length-1 is the ']' character
        my $nested_depth = 1;

        while ($nested_depth && ($index > -1) ) {

            my $char = substr $description, $index, 1;
            CHAR: {
                ($char eq '[') && do {
                    $nested_depth--;
                    last CHAR;
                };
                ($char eq ']') && do {
                    $nested_depth++;
                    last CHAR;
                };
            }
            $index--;

        }

        croak qq{unexpected header:\n$header\n} if ($nested_depth);

        # At this point, $index is one to the left of the closing '[' (or opening, if reading L->R)

        $taxon_string   = substr $description, $index + 2, length($description) - 3 - $index;
        $desc_remainder = substr $description, 0, $index;

=pod
        example $description string:

       'hypothetical protein Celgi_0332 [[Cellvibrio] gilvus ATCC 13127]'
        0123456789012345678901234567890123456789012345678901234567890123
        0         1         2         3         4         5         6

        string length is 64

        end of taxon string is 64 - 2 = 62; counting down from 62,
        $nested_depth is incremented to 2 at position 44; decremented to
        1 at 33; and to 0 at 32; so the start position of the taxon
        string is 32 + 1 = 33, and the length is 1 + 62 - 33 = 30.
=cut

    }

    return ($seq_id, $desc_remainder, $taxon_string);
}

sub output_fasta_record {

    my $fasta_record_lref = shift;
    my $fh                = shift || \*STDOUT;

    my @fasta_record = @{$fasta_record_lref};

    print $fh @fasta_record;
}

