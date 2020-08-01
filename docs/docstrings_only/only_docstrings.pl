#  generates .py files from .py(x)
#  source code with docstring only

my @files = glob('../fathon/*.py*');
foreach my $file (@files)
{
    open my $in,'<',$file or die "Can't open $file\n";
    
    my $file_name = (split('/', $file))[-1];
    my $file_name_no_ext = (split('\.', $file_name))[0];
    my $out_file = "./docstrings_only/fathon/$file_name_no_ext.py";
    open my $out,'>',$out_file or die "Can't open $out_file\n";

    $write_to_file = 0;
    $class_docs = 0;
    while(<$in>)
    {
        chomp($_);
        
        if($file_name_no_ext =~ m/__init__/)
        {
            next if m/^__/;
            if($_ =~ m/\"{3}.+/)
            {
                print $out "from .fathonUtils import *\n";
            }
            print $out $_."\n";
        }
        else
        {
            if($_ =~ m/^(from|import)/)
            {
                print $out $_."\n";
            }
            
            if($_ =~ m/^(\s+)?(def|cpdef)/)
            {
                $class_docs = 0;
                
                if($_ =~ m/__init__/)
                {
                    print $out $_."\n";
                    $_ =~ s/def.*/\tpass/g;
                    print $out $_."\n";
                    print $out "\n";
                }
                else
                {
                    $_ =~ s/cpdef/def/g;
                    $_ =~ s/, .*\[.+\]/,/g;
                    $_ =~ s/, [^,]+ /, /g;
                    print $out $_."\n";
                }
            }
            
            if($_ =~ m/cdef class/)
            {
                $class_docs = 1;
                
                $_ =~ s/cdef //g;
                print $out $_."\n";
            }
            
            if($_ =~ m/^\#(?=\s+)/)
            {
                print $out $_."\n";
            }
            
            if($_ =~ m/\"{3}.+/)
            {
                $write_to_file = 1;
            }
            elsif($_ =~ m/\"{3}$/)
            {
                print $out $_."\n";
                if($class_docs == 0)
                {
                    $_ =~ s/\"{3}/return 0/g;
                    print $out $_."\n";
                }
                print $out "\n";
                $write_to_file = 0;
            }
            if($write_to_file == 1)
            {
                print $out $_."\n";
            }
        }
    }

    close $in;
    close $out;
    
}
