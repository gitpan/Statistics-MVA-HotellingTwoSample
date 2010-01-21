package Statistics::MVA::HotellingTwoSample;

use base Statistics::MVA;
use warnings;
use strict;
use Carp;
use Statistics::Distributions qw/fprob/;

=head1 NAME

Statistics::MVA::HotellingTwoSample - Two-Sample Hotelling's T-Square Test Statistic.

=cut
=head1 VERSION

This document describes Statistics::MVA::HotellingTwoSample version 0.0.1

=cut
=head1 SYNOPSIS

    use Statistics::MVA::HotellingTwoSample;

    # we have two groups of data each with 3 variables and 10 observations - example data from http://www.stat.psu.edu/online/courses/stat505/data/insect.txt
    my $data_X = [
        [qw/ 191 131 53/],
        [qw/ 185 134 50/],
        [qw/ 200 137 52/],
        [qw/ 173 127 50/],
        [qw/ 171 128 49/],
        [qw/ 160 118 47/],
        [qw/ 188 134 54/],
        [qw/ 186 129 51/],
        [qw/ 174 131 52/],
        [qw/ 163 115 47/],
    ];

    my $data_Y = [
        [qw/ 186 107 49/],
        [qw/ 211 122 49/],
        [qw/ 201 144 47/],
        [qw/ 242 131 54/],
        [qw/ 184 108 43/],
        [qw/ 211 118 51/],
        [qw/ 217 122 49/],
        [qw/ 223 127 51/],
        [qw/ 208 125 50/],
        [qw/ 199 124 46/],
    ];
    
    # Create a Statistics::MVA::HotellingTwoSample object and pass the data as two Lists-of-Lists within an anonymous array.
    my $mva = Statistics::MVA::HotellingTwoSample->new([ $data_X, $data_Y ], {standardise => 0});

    # Generate results and print a report to STDOUT by calling hotelling_two_sample in void context.
    $mva->hotelling_two_sample;

    # Call hotelling_two_sample in LIST-context to access the results directly.
    my ($T2, $F, $pval, $df1, $df2) = $mva->hotelling_two_sample;

=cut
=head1 DESCRIPTION

Hotelling's T-square statistics is a generalization of Student's t statistic that is used for multivariate hypothesis testing.

=cut

#/ for the matrixreal objects of CV mats: $self->[0][$i][0];

#/ for a cv matrix $self->[0][$i][0][0];

#/ for raw data: self->[3][$i]

use version; our $VERSION = qv('0.0.1');

#/ psu prints covariance matrices too...

sub hotelling_two_sample {

    my $self = shift;

    my $k = $self->[1];
    croak qq{\nThis is a two-sample test - you must give two samples} if $k != 2;  

    my $n_x = $self->[0][0][1];
    my $n_y = $self->[0][1][1];

    #y just variable number - again will need to check its equal for sample 2 - Statistics::MVA already checked p is same for all matrices...
    my $p = $self->[2];

    #y Just averages of for each variable
    my @bar_x = &_bar($self->[3][0]);
    my @bar_y = &_bar($self->[3][1]);

    #y covariance matrices - this is already done!!! as part of MVA object creation
    my $V_x = $self->[0][0][0];
    my $V_y = $self->[0][1][0];
        
    #/ pooled is the matrix addition sum of the individual cv_matrices each multiplied by n_i-1 divided by (n_a + n_b - 2): ((99*Vx)+(Vy*99))/(198)
    #b $Vp <- (($n - 1) * $V + ($ny - 1) * $Vy) / ($n + $ny - 2)
    
    my $V_x_adj = $V_x->shadow;
    my $V_y_adj = $V_y->shadow;
    my $V_p = $V_x->shadow; 

    $V_x_adj->multiply_scalar($V_x,$n_x-1);
    $V_y_adj->multiply_scalar($V_y,$n_y-1);
    $V_p->add($V_x_adj,$V_y_adj); 

    my $divisor = 1 / ($n_x + $n_y - 2);
    $V_p->multiply_scalar($V_p,$divisor);
    
    # subtract means...
    my @bar_diff = map { $bar_x[$_] - $bar_y[$_] } (0..$#bar_x);

    my $bar_diff_mat = Math::MatrixReal->new_from_cols([[@bar_diff]]);

    my $dim;
    my $x;
    my $B_matreal;

    #/ using matrixreal - SOLVE: generic function solves the equation 'a %*% x = b' for 'x', where 'b' can be either a vector or a matrix.
    my $LR = $V_p->decompose_LR();
    if (!(($dim,$x,$B_matreal) = $LR->solve_LR($bar_diff_mat))) { croak qq{\nsystem has no solution} }
    #print qq{\n\nhere is the matrixreal solution\n}, $x;

    #/ using Cephes - thus perhaps put all in Cephes - i.e. just addition/subtraction and scalar multiplication until now..
    #my $M = Math::Cephes::Matrix->new($V_p->[0]);
    #my $B = [@bar_diff];
    #my $solution = $M->simq($B);
    #print qq{\nhere is the cephes solution @{$solution}\n};
    
    #  ((t(dif in means) %*% inverse(Vp, dif.means)) * ( (100+100-6-1) * ((100*100)/(100+100) ) * (198) * 6 ) / ((100+100-2) * 6 * 193))
    #y crossprod * ( (100+100-6-1) * ((100*100)/(100+100) ) * (198) * 6 ) / ((100+100-2) * 6 * 193))
    my $crossprod = ~$bar_diff_mat * $x;

    my $df1 = $p;
    my $df2 = $n_x + $n_y-$p-1;
    my $other = $n_x + $n_y - 2;

    # F = (((100+100-6-1)*2412.45))/(198*6) 
    #   = ((n_a + n_b - k - 1)*T^2)/((n_a + n_b - 2) * k)
    # T^2 = (391.9216*(100+100-2)*6)/(100+100-6-1)
    #     = ( F * k * (n_a + n_b - 2) ) / (n_a + n_b - 6 - 1)

    my $T2 = ($crossprod->[0][0][0] * $df2 * (($n_x * $n_y)/($n_x+$n_y)) * $other * $p) / (($n_x+$n_y-2) * $p * $df2);
    #my $T2 = ($crossprod->[0][0][0] * ($n_x + $n_y-$p-1) * (($n_x * $n_y)/($n_x+$n_y)) * ($n_x + $n_y - 2) * $p) / (($n_x+$n_y-2) * $p * ($n_x + $n_y-$p-1));
    
    my $F = ( $df2 * $T2 ) /  ($other * $df1);
    my $pval = &fprob($df1,$df2,$F);

    $pval = $pval < 1e-8 ? q{< 1e-8} : $pval;

    # t2=(ybar1-ybar2)`*inv(spool*(1/n1+1/n2))*(ybar1-ybar2);
    # f=(n1+n2-k-1)*t2/k/(n1+n2-2);
    
    if ( !wantarray ) { print qq{\nT^2 = $T2\nF = $F \ndf1 = $df1\ndf2 = $df2 \np.value = $pval}; }
    else { return ( $T2, $F, $pval, $df1, $df2 ) }
    
    return;
}

sub _bar {

    my $lol = shift;
    my $rows = scalar @{$lol};
    my $cols = scalar @{$lol->[0]};
    my @bar;
    # not necessary
    #$#bar = scalar @{$lol->[0]}-1;

    for my $col (0..$cols-1) {
        my $sum = 0;
        for my $row (0..$rows-1) {
            $sum += $lol->[$row][$col];
        }
        push @bar, $sum/$rows;
    }

    $, = q{, };

    return @bar;
}

1; # Magic true value required at end of module

__END__

=head1 DEPENDENCIES

'Statistics::MVA' => '0.0.1',
'Carp' => '1.08',
'Statistics::Distributions' => '1.02',

=cut
=head1 BUGS AND LIMITATIONS

Let me know.

=cut
=head1 AUTHOR

Daniel S. T. Hughes  C<< <dsth@cantab.net> >>

=cut
=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010, Daniel S. T. Hughes C<< <dsth@cantab.net> >>. All rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=cut
=head1 DISCLAIMER OF WARRANTY

Because this software is licensed free of charge, there is no warranty
for the software, to the extent permitted by applicable law. Except when
otherwise stated in writing the copyright holders and/or other parties
provide the software "as is" without warranty of any kind, either
expressed or implied, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose. The
entire risk as to the quality and performance of the software is with
you. Should the software prove defective, you assume the cost of all
necessary servicing, repair, or correction.

In no event unless required by applicable law or agreed to in writing
will any copyright holder, or any other party who may modify and/or
redistribute the software as permitted by the above licence, be
liable to you for damages, including any general, special, incidental,
or consequential damages arising out of the use or inability to use
the software (including but not limited to loss of data or data being
rendered inaccurate or losses sustained by you or third parties or a
failure of the software to operate with any other software), even if
such holder or other party has been advised of the possibility of
such damages.
=cut
