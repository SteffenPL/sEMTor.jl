function write_sheet!(df, xf, prop, type;
                            cols=:time,
                            col_name=string(cols),
                            prop_name=string(prop),
                            df_mean=missing,
                            max_reps=500)

    if isempty(df[!,prop])
        return nothing
    end

    sheetname = prop_name * "_" * type
    if length(sheetname) > 31
        sheetname = replace(sheetname, "apical" => "a")
        sheetname = replace(sheetname, "basal" => "b")
        sheetname = replace(sheetname, "nuclei" => "n")
    end

   if "Sheet1" == XLSX.sheetnames(xf)[1]
       sheet = xf[1]
       XLSX.rename!(sheet, sheetname)
   else
       sheet = XLSX.addsheet!(xf, sheetname)
   end


   rep_offset = ismissing(df_mean) ? 0 : 2

   sheet[1,2 + rep_offset] = "rep"
   sheet[2,1] = col_name

   sheet["A3", dim=1] = unique(df[!,cols])
   col = 2
   dff = filter( row -> row.rep <= max_reps, df)
   for (key2, dfrep) in pairs(groupby(dff, :rep))
       sheet[2,col + rep_offset] = key2.rep
       sheet[3,col + rep_offset, dim=1] = collect(dfrep[!,prop])
       col += 1
   end

   if !ismissing(df_mean)

        sheet[2,2] = "mean"
        sheet[3,2, dim=1] = collect(df_mean[!,Symbol(prop,:_mean)])
        sheet[2,3] = "std"
        sheet[3,3, dim=1] = collect(df_mean[!,Symbol(prop,:_std)])
   end

   nothing
end


function export_statistics_to_xlsx(stats, stats_div, stats_mean, fn; max_reps = 500)
    if !occursin(".xlsx", fn)
        fn = fn * ".xlsx"
    end

    XLSX.openxlsx(fn, mode="w") do xf

        for (cell_type, df) in stats
            cols = filter( x -> !(x in ["rep", "time", "label", "type"]), names(df))
            for col in cols
                write_sheet!(df, xf, col, cell_type, df_mean = stats_mean[cell_type]; max_reps)
            end
        end


        for (cell_type, df_div) in stats_div
            cols = filter( x -> !(x in ["rep", "time", "label", "type"]), names(df_div))
            for col in cols
                write_sheet!(df_div, xf, col, cell_type,
                    cols=:label, col_name="cell label", prop_name="div_" * string(col); max_reps)
            end
        end
    end
    nothing
end
