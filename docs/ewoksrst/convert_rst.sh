for file in *.rst; do
    pandoc "$file" -f rst -t markdown -o "../ewoks_specs/${file%.rst}.md"
done

