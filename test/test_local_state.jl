using ITensors

@testset "local_state" begin

    @testset "local_state_index / local_state_string are inverses" begin
        for (sitetype, names) in [
            ("S=1/2",    ["Up", "Dn"]),
            ("tJ",       ["Emp", "Up", "Dn"]),
            ("Electron", ["Emp", "Up", "Dn", "UpDn"]),
        ]
            s = ITensors.siteind(sitetype)

            @testset "$(sitetype): string → index → string" begin
                for name in names
                    @test local_state_string(s, local_state_index(s, name)) == name
                end
            end

            @testset "$(sitetype): index → string → index" begin
                for n in 1:length(names)
                    @test local_state_index(s, local_state_string(s, n)) == n
                end
            end
        end
    end

end
