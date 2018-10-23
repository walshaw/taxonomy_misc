#!/usr/bin/perl

use strict;
use warnings;

use Carp;

use Getopt::Long;
use List::Util qw( first );

# These are set (optionally in some cases) by the user
my $taxonomy_nodes_file;
my $max_iterations   = 2;
my $report_frequency = 10000;
my $sorted_outfile;
#(end)

my $delimiter        = q{|};
#my $delimiter_qr     = qr{$delimiter}xms;
my $delimiter_quoted = quotemeta $delimiter;

#print qq{$delimiter_qr\n}; exit;

my @hierarchy = (
        'root'        ,
	'superkingdom', 'kingdom'       , 'subkingdom'      ,
	'superphylum' , 'phylum'        , 'subphylum'       ,
	'superclass'  , 'class'         , 'subclass'        , 'infraclass',
	'superorder'  , 'order'         , 'suborder'        , 'infraorder', 'parvorder',
	'superfamily' , 'family'        , 'subfamily'       , 
                        'tribe'         , 'subtribe'        ,
                        'genus'         , 'subgenus'        ,
                        'species group' , 'species subgroup',
                        'species'       , 'subspecies'      ,
                        'varietas',
                        'forma',
);

my %ordinal_rank_of;
for my $i ( 0 .. $#hierarchy ) { $ordinal_rank_of{$hierarchy[$i]} = $i; }
$ordinal_rank_of{'no rank'} = -1;

my $usage = qq{Usage:\n\n$0 [ -nodesfile ] TAXONOMY_NODES_FILE [ options ]

Options:
	-iterations	integer		MAX_ITERATIONS     (default $max_iterations)
	-reportfrq	integer		REPORT_FREQUENCY   (default $report_frequency)
	-sortednodes	string		SORTED_OUTPUT_FILE (default - none)

Probably want to use nodes.dmp from NCBI Taxonomy, as TAXONOMY_NODES_FILE.

The file is repeatedly read* until no further nodes are added to the clade,
or until the file has been read MAX_ITERATIONS times.

The point is that it cannot be assumed that the nodes are listed in a
'hierarchical order', i.e. it cannot be assumed that every node is listed
after its parent, grandparent and all ancestral nodes.

*Actually, although the file may be large - currently nodes.dmp is 1.3 million
lines long - it is still read only once (only the first two fields are read,
which both contain integer values), and then the data iterated over.
};

croak qq{$usage} if !(GetOptions(
        "nodesfile=s"   => \$taxonomy_nodes_file,
        "iterations=i"  => \$max_iterations,
        "reportfrq=i"   => \$report_frequency,
        "sortednodes=s" => \$sorted_outfile,
    )
);

$taxonomy_nodes_file ||= shift;

croak $usage if (!defined $taxonomy_nodes_file);

croak qq{file '$sorted_outfile' already exists; please delete it first if you are sure you want to overwrite it}
    if ($sorted_outfile && -f $sorted_outfile);

my $universal_clade_node = 1; # the ancestor of everything else; in the nodes.dmp
                              # file, the taxonomic rank is 'no rank' (N.B. many
                              # other, much more specific, nodes also have this
                              # rank) and its parent is *itself*

my @parent_of;     # indices correspond to node IDs; values are node IDs of parent
my @rank_of;       # indices correspond to node IDs; values are strings (each is the name of a rank)
my @rank_order_of; # indices correspond to node IDs; values are integers, each representing a rank
                   # (see %ordinal_rank_of)
my @remainder_of;  # indices correspond to node IDs; values are strings

my %rank_freq; # keys strings (rank names); values are integers (incidence in file)

=pod

    A comprehensive list of the taxonomic ranks used in the nodes.dmp file can be easily
    obtained:
    $ perl -ne '@fields = split /\s+\|\s+/; print "$fields[2]\n"; ' nodes.dmp | sort -u > ranks.lst 

    
=cut

my %parent_rank_of_rank; # keys are strings (taxonomic levels); values are refs to hashes,
                         # whose keys are the parent taxonomic levels, and whose values are
                         # integers (number of times the specified parent is the parent rank
                         # of the first rank


open my $fh, '<', $taxonomy_nodes_file or croak qq{couldn't read nodes file '$taxonomy_nodes_file'};

my $n_lines = 0;

my $n_uncharacterized_parents = 0; # basically, some nodes are listed before their parent node;
                                   # so the characteristics of the parent will be unknown at 
                                   # the time the child node is read; by having 2 iterations,
                                   # this number can be reduced to 0


LINE:
while (defined (my $line = <$fh>)) {

    my @fields = split m{ \s* $delimiter_quoted \s* }xms, $line, 4;

    my ($node_id, $parent_id, $rank, $remainder) = @fields[0,1,2,3];

    # $remainder should begin with a 
    # tab or a non-whitespace field; the tab will have been lost by the split if
    # the field was empty;
    # i.e. an empty field should be: '|		|' not '|	|'

    $remainder =~ s{ \A ($delimiter_quoted) }{\t$1}xms;
    $parent_of[$node_id]    = $parent_id;
    $rank_of[$node_id]      = $rank;
    $remainder_of[$node_id] = $remainder; # the rest of the line following the 3 fields of interest

    if (defined $rank_of[$parent_id]) {
        $parent_rank_of_rank{$rank}{$rank_of[$parent_id]}++;
    }
    else {
        #print STDERR qq{node $node_id has parent $parent_id; parent has no known rank yet\n};
        $n_uncharacterized_parents++;
    }

    $rank_freq{$rank}++;

#    if (first { $_ == $parent_id } @clade_nodes) {
#        push @clade_nodes, $node_id;
        ###print STDERR qq{iteration 1: adding $node_id (parent is $parent_id)\n};
#    }

    if ( !(++$n_lines % $report_frequency) ) {
            
        print STDERR qq{$n_lines lines read\n};
    }

}

print STDERR qq{$n_lines lines read\n};

close $fh or croak qq{couldn't close nodes file '$taxonomy_nodes_file'};

print STDERR qq{completed first iteration; $n_uncharacterized_parents },
             qq{nodes have parents that were uncharacterized (listed afterwards)\n};

my $n_iterations = 1;

ITERATION:
while ($n_iterations < $max_iterations) {

    $n_iterations++;

    my $n_uncharacterized_parents = 0;

    NODE:
    for my $node_id (0 .. $#parent_of) {

        next NODE if (!defined $parent_of[$node_id]); # some node_ids are unused in the nodes.dmp file

        my $parent_id = $parent_of[$node_id];
        my $rank      = $rank_of[$node_id];

        if (defined $rank_of[$parent_id]) {
            $parent_rank_of_rank{$rank}{$rank_of[$parent_id]}++;
        }
    else {
        print STDERR qq{node $node_id has parent $parent_id; parent has no known rank yet\n};
        $n_uncharacterized_parents++;
        }
    }

    print STDERR qq{completed iteration $n_iterations; $n_uncharacterized_parents },
             qq{nodes have parents that were uncharacterized (listed afterwards)\n};

}

print qq{Overall frequencies of each rank:\n\n};

my @sorted_ranks = sort { $rank_freq{$a} <=> $rank_freq{$b} } (keys %rank_freq);

for my $rank (@sorted_ranks) {

    print qq{\t$rank\t$rank_freq{$rank}\n};

}

print qq{\n\n};

for my $rank (sort keys %parent_rank_of_rank) {

    print qq{rank '$rank' has parent ranks with these frequencies:\n};

    for my $parent_rank (sort keys %{$parent_rank_of_rank{$rank}}) {

        print qq{\t$parent_rank\t$parent_rank_of_rank{$rank}{$parent_rank}\n};

    }
    print qq{\n};

}

if ($sorted_outfile) {

    

    print STDERR qq{\n\ncreating output file containing sorted (by descending rank) nodes\n\n};

    open my $fh, '>', $sorted_outfile or croak qq{couldn't create output file $sorted_outfile};

    print STDERR qq{assigning ordinal rank to all nodes\n\n};

    my @used_ids;

    NODE:
    for my $node_id (0 .. $#parent_of) {

        next NODE if (!defined $parent_of[$node_id]); # some node_ids are unused in the nodes.dmp file

        # Because a lot of node IDs are unused in the NCBI listing (probably old records expunged)
        # a lot of elements of @parent_of are therefore undef. That makes them less than suitable
        # for input to a sort; even though the undef values can be handled conditionally, this
        # would be a bit wasteful in terms of time - it is a very big list after all. So a new
        # list is made, of just the used numbers.

        push @used_ids, $node_id;


        my $rank = $rank_of[$node_id];
        croak qq{unknown rank:\t'$rank' (node $node_id)} if (!defined $ordinal_rank_of{$rank});

        my $rank_order = $ordinal_rank_of{$rank}; # if $rank eq 'no rank' then this is -1

        my $test_rank_no  = $rank_order;
        my $rank_offset = 0;
        my $next_node   = $node_id;
        ###my $rank_name;

        if ($node_id == $universal_clade_node) {
            $rank_order_of[$node_id] = 0; # universal ancestral, i.e. 'root' node; this will
                                          # also be 'no rank' and so by default it would have
                                          # an ordinal rank of -1
            next NODE;
        }

        while (($test_rank_no < 0) && ($next_node != 1)) { # it's a 'no rank'; so see what its parent is; and if necessary,
                                  # its grandparent etc (parent could also be 'no rank')
            $next_node    = $parent_of[$next_node];
            my $next_rank = $rank_of[$next_node];
            $test_rank_no = $ordinal_rank_of{$next_rank};
            $rank_offset++;
            print STDERR qq{$node_id ('no rank'): testing next ancestral node $next_node ($next_rank) ($rank_offset level(s) above)\n};
        }

        my $inferred_rank_no = $test_rank_no + $rank_offset;

        if ($rank_offset) { # only report the 'no rank' cases
            print STDERR qq{assigning (notional) rank order $inferred_rank_no },
                         qq{(notionally equivalent to $hierarchy[$inferred_rank_no]) },
                         qq{to node $node_id ($rank)\n\n};
        }

        $rank_order_of[$node_id] = $inferred_rank_no; # unless the node's original rank is 'no rank' then this
                                                      # value is identical to the ordinal of the original rank
    }

    print STDERR qq{highest rank no. is $#parent_of; total number of unique used ranks is }, scalar @used_ids, qq{\n\n};

    print STDERR qq{sorting node IDs in ascending order of rank (or notional rank)\n\n};

    my @sorted_used_ids = sort { $rank_order_of[$a] <=> $rank_order_of[$b] } @used_ids;

    for my $node_id (@sorted_used_ids) {

        # $remainder_of[] element will end in a \n character
        print $fh qq{$node_id\t|\t$parent_of[$node_id]\t|\t$rank_of[$node_id]\t|\t$remainder_of[$node_id]};

    }

    close $fh or croak qq{couldn't close output file $sorted_outfile};

    print STDERR qq{sorting complete.\n\n};
}


