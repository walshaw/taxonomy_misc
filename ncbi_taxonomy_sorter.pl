#!/usr/bin/perl

use strict;
use warnings;

use Carp;

use List::Util qw( first );

my $default_max_iterations = 100;

my $delimiter        = q{|};
#my $delimiter_qr     = qr{$delimiter}xms;
my $delimiter_quoted = quotemeta $delimiter;

#print qq{$delimiter_qr\n}; exit;

my $usage = qq{Usage:\n\n$0 TAXONOMY_NODES_FILE ANCESTRAL_NODE_ID [ MAX_ITERATIONS  [ REPORT_FREQUENCY ] ]

Probably want to use nodes.dmp from NCBI Taxonomy, as the first arg.
ANCESTRAL_NODE_ID is the node at the very top of the tree of interest.
The file is repeatedly read* until no further nodes are added to the clade,
or until the file has been read MAX_ITERATIONS times.

The point is that it cannot be assumed that the nodes are listed in a
'hierarchical order', i.e. it cannot be assumed that every node is listed
after its parent, grandparent and all ancestral nodes.

*Actually, although the file may be large - currently nodes.dmp is 1.3 million
lines long - it is still read only once (only the first two fields are read,
which both contain integer values), and then the data iterated over.
};

my $taxonomy_nodes_file = shift or croak $usage;

my $ancestral_node      = shift or croak $usage;

my $max_iterations      = shift || $default_max_iterations; # so this can't be 0

my $report_frequency    = shift || 10000;

my $universal_clade_node = 1; # the ancestor of everything else; in the nodes.dmp
                              # file, the taxonomic rank is 'no rank' (N.B. many
                              # other, much more specific, nodes also have this
                              # rank) and its parent is *itself*

my @parent_of; # indices correspond to node IDs; values are node IDs of parent
my @rank_of;   # indices correspond to node IDs; values are strings

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

my @clade_nodes = ($ancestral_node);

my $n_lines = 0;

LINE:
while (defined (my $line = <$fh>)) {

    my @fields = split m{ \s* $delimiter_quoted \s* }xms, $line;

    my ($node_id, $parent_id, $rank) = @fields[0,1,2];

    $parent_of[$node_id] = $parent_id;
    $rank_of[$node_id]   = $rank;
    if (defined $rank_of[$parent_id]) {
        $parent_rank_of_rank{$rank}{$rank_of[$parent_id]}++;
    }

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

for my $rank (sort keys %parent_rank_of_rank) {

    print qq{rank '$rank' has parent ranks with these frequencies:\n};

    for my $parent_rank (sort keys %{$parent_rank_of_rank{$rank}}) {

        print qq{\t$parent_rank\t$parent_rank_of_rank{$rank}{$parent_rank}\n};

    }
    print qq{\n};

}


