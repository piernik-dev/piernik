for i in jenkins/gold_configs/*.config ; do
	case $(basename $i) in
        (mcrtest_CRESP.config|mcrtest_new.config|mcrwind.config|MHDsedovAMR.config|resist.config|streaming_instability.config) ;; # te są zaimplementowane
        (*) ./jenkins/gold_test.sh $i ;;
    esac
done

FAKE_DIR=jenkins/workspace/no_tests/
[ -e $FAKE_DIR ] && rm -rf $FAKE_DIR
