filename = 'Dataset_test_PC';
sheet = 'Quarterly';
range = 'B2:DA300';
do_truncation = 1; %Do truncate data.
[dataset, var_names] = read_data2(filename, sheet, range, do_truncation);
dataset = real(dataset);
date_start = dataset(1,1);
dataset = dataset(:,2:end); %Removing time befor PC analysis

pc = get_principal_components(dataset);

PC_BIG = nan(size(dataset,1),size(datset,2))

hold on
for i=1:5
      plot(pc(:,i))
end