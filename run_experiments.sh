#!/bin/bash

function prepare_options {
  cat <<EOF > options.$1.txt;
--gen-method $1
--num-trials $2
--n-max $3
--n-min $3
--verbosity $4
--print-path $5
--print-tree $5
--measure-mle $6
--test-confidence $7
--lambda 1
--k $8
--p $9
--dump-numbers ${10}
--use-random-weights false
--num-threads 1
EOF
}

function run_test {
  num_args=$#;

  if [ $num_args -lt 5 ]; then
    echo "Incorrect number of arguments: " $num_args;
  fi

  distribution=$1;
  num_trials=$2;
  num_nodes=$3;
  verbose=$4;
  print_tree=$5;

  if [ $num_args -gt 5 ]; then
    measure_mle=$6;
    test_confidence=$7;
  else
    measure_mle=true;
    test_confidence=false;
  fi

  if [ $num_args -gt 7 ]; then
    k=$8;
    p=$9;
  else
    k=2;
    p=0.5;
  fi

  if [ $num_args -gt 9 ]; then
    dump_numbers=${10};
  else
    dump_numbers=false;
  fi

  prepare_options ${distribution} \
                  $num_trials \
                  $num_nodes \
                  $verbose \
                  $print_tree \
                  ${measure_mle} \
                  ${test_confidence} \
                  ${k} \
                  ${p} \
                  ${dump_numbers}
  result_file=result-$distribution-$num_trials-$num_nodes-$k.txt;
  ./harness @options.$distribution.txt > $result_file;
}

for arg in "$@"; do
  case ${arg} in
    --toy-Binomial)
    run_test Binomial 1 11 2 true true false;
    ;;

    --toy-Geometric)
    run_test Geometric 1 11 2 true true false;
    ;;

    --toy-Poisson)
    run_test Poisson 1 11 2 true true false;
    ;;

    --toy-Binary-0-2)
    run_test Binary-0-2 1 11 2 true true false;
    ;;

    --toy-Binary-0-1-2)
    run_test Binary-0-1-2 1 11 2 true true false;
    ;;

    --toy-all)
    for method in Binomial Poisson Geometric Binary-0-2 Binary-0-1-2; do
      run_test $method 1 11 2 true true false;
    done
    ;;

    --test-Binomial)
    for trials in 1000 5000 10000; do
      for n in 1001 2001 4001; do
        run_test Binomial $trials $n 1 false true false;
      done
    done
    ;;

    --test-Geometric)
    for trials in 1000 5000 10000; do
      for n in 1001 2001 4001; do
        run_test Geometric $trials $n 1 false true false;
      done
    done
    ;;

    --test-Poisson)
    for trials in 1000 5000 10000; do
      for n in 1001 2001 4001; do
        run_test Poisson $trials $n 1 false true false;
      done
    done
    ;;

    --test-Binary-0-2)
    for trials in 1000 5000 10000; do
      for n in 1001 2001 4001; do
        run_test Binary-0-2 $trials $n 1 false true false;
      done
    done
    ;;

    --test-Binary-0-1-2)
    for trials in 1000 5000 10000; do
      for n in 1001 2001 4001; do
        run_test Binary-0-1-2 $trials $n 1 false true false;
      done
    done
    ;;

    --test-all)
    for method in Binomial Poisson Geometric Binary-0-2 Binary-0-1-2; do
      for trials in 1000 5000 10000; do
        for n in 1001 2001 4001; do
          run_test $method $trials $n 1 false true false;
        done
      done
    done
    ;;

    --test-confidence)
    for method in Binomial Poisson Geometric Binary-0-2 Binary-0-1-2; do
      for trials in 5000 10000; do
        if [ "$method" == 'Binomial' ]; then
          echo "run_test $method $trials 1001 0 false false true 2 0.5;"
          run_test $method $trials 1001 0 false false true 2 0.5;
          #echo "run_test $method $trials 1001 0 false false true 4 0.25;"
          #run_test $method $trials 1001 0 false false true 4 0.25;
          #echo "run_test $method $trials 1001 0 false false true 8 0.125;"
          #run_test $method $trials 1001 0 false false true 8 0.125;
        else
          echo "run_test $method $trials 1001 0 false false true;"
          run_test $method $trials 1001 0 false false true;
        fi
      done
    done
    ;;

    --dump-numbers-Binomial)
    for method in Binomial; do
      filename=${method}-big-file.txt
      for i in `seq 1 1000`; do 
        rm -rf result-${method}-*.txt;
        echo "${method} ${i}";
        run_test $method 1000 1001 0 false false false 2 0.5 true;
        cat result-${method}-* >> ${filename}
      done
    done
    ;;

    --dump-numbers-Poisson)
    for method in Poisson; do
      filename=${method}-big-file.txt
      for i in `seq 1 1000`; do 
        rm -rf result-${method}-*.txt;
        echo "${method} ${i}";
        run_test $method 1000 1001 0 false false false 2 0.5 true;
        cat result-${method}-* >> ${filename}
      done
    done
    ;;

    --dump-numbers-Geometric)
    for method in Geometric; do
      filename=${method}-big-file.txt
      for i in `seq 1 1000`; do 
        rm -rf result-${method}-*.txt;
        echo "${method} ${i}";
        run_test $method 1000 1001 0 false false false 2 0.5 true;
        cat result-${method}-* >> ${filename}
      done
    done

    ;;

    --dump-numbers-Binary-0-2)
    for method in Binary-0-2; do
      filename=${method}-big-file.txt
      for i in `seq 1 1000`; do 
        rm -rf result-${method}-*.txt;
        echo "${method} ${i}";
        run_test $method 1000 1001 0 false false false 2 0.5 true;
        cat result-${method}-* >> ${filename}
      done
    done
    ;;

    --dump-numbers-Binary-0-1-2)
    for method in Binary-0-1-2; do
      filename=${method}-big-file.txt
      for i in `seq 1 1000`; do 
        rm -rf result-${method}-*.txt;
        echo "${method} ${i}";
        run_test $method 1000 1001 0 false false false 2 0.5 true;
        cat result-${method}-* >> ${filename}
      done
    done
    ;;

    *)
    echo "run_experiments: options are"
echo "--test-all: test all the distributions (Poisson, Binomial, Geometric)"
echo "--test-Binomial: test Binomial distribution"
echo "--test-Poisson: test Poisson distribution"
echo "--test-Geometric: test Geometric distribution"
echo "--test-Binary-0-2: test Binary-0-2 distribution"
echo "--test-Binary-0-1-2: test Binary-0-1-2 distribution"
echo "--test-Geometric: test Geometric distribution"
echo "--toy-all: toy all the distributions (Poisson, Binomial, Geometric)"
echo "--toy-Binomial: toy Binomial distribution"
echo "--toy-Poisson: toy Poisson distribution"
echo "--toy-Geometric: toy Geometric distribution"
echo "--toy-Binary-0-2: toy Binary-0-2 distribution"
echo "--toy-Binary-0-1-2: toy Binary-0-1-2 distribution"
echo "--test-confidence: testing Confidence intervals"
    ;;
  esac
done
